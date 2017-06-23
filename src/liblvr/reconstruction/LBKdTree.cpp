#include <stdio.h>

#include <iostream>
#include "lvr/reconstruction/LBKdTree.hpp"


// Static variables
ctpl::thread_pool* LBKdTree::pool = new ctpl::thread_pool(8);
int LBKdTree::st_num_threads = 8;
int LBKdTree::st_depth_threads = 3;

/// Public

LBKdTree::LBKdTree( LBPointArray& vertices, int num_threads) {
    st_num_threads = num_threads;
    st_depth_threads = static_cast<int>(log2(st_num_threads));
    pool = new ctpl::thread_pool(st_num_threads);
    this->generateKdTree(vertices);
}

void LBKdTree::generateKdTree(LBPointArray &vertices) {
    struct LBPointArray indices_sorted[ vertices.dim ];
    struct LBPointArray values_sorted[ vertices.dim ];

    // 3 dims in threads!
    //this->p.push(addFloatPointerRec, vec_a.get(), vec_b.get(), vec_res_thr.get(), 0, vec_size);

    for(int i=0; i< vertices.dim; i++)
    {
        pool->push(generateAndSort, vertices, &(indices_sorted[i]), &(values_sorted[i]), i);
    }
    
    pool->stop(true);
    pool = new ctpl::thread_pool(st_num_threads);
    
    this->generateKdTreeArray(vertices, indices_sorted, vertices.dim, this->kd_tree);

    for(int i=0; i<vertices.dim;i++)
    {
        //free(indices_sorted[i].elements);
        free(values_sorted[i].elements);
    }
}

LBPointArray LBKdTree::getKdTreeArray() {
    return this->kd_tree;
}

/// Private

void LBKdTree::generateKdTreeArray(LBPointArray& V, LBPointArray* sorted_indices, int max_dim, LBPointArray& kd_tree) {

    int size;
    int max_tree_depth;

    max_tree_depth = static_cast<int>( log2f(V.width - 1 ) + 2.0 ) ;

    if (V.width == 1)
    {
        max_tree_depth = 1;
    }

    size = V.width * 2 - 1;

    generatePointArray(kd_tree, size, 1);

    //start real generate
    generateKdTreeRecursive(0, V, sorted_indices, 0, max_dim, kd_tree, size, max_tree_depth, 0, 0);

    //pool->push(generateKdTreeRecursive, V, sorted_indices, 0, max_dim, kd_tree, size, max_tree_depth, 0, 0);

    pool->stop(true);
    pool = new ctpl::thread_pool(st_num_threads);
}

void LBKdTree::generateKdTreeRecursive(int id, LBPointArray& V, LBPointArray sorted_indices[], int current_dim, int max_dim, LBPointArray& kd_tree, int size, int max_tree_depth, int position, int current_depth) {
        
    int left = position*2+1;
    int right = position*2+2;


    if( right > size-1 || left > size-1 )
    {

        kd_tree.elements[position] = sorted_indices[current_dim].elements[0];

    } else {
        /// split sorted_indices
        int indices_size = sorted_indices[current_dim].width;

        int v = pow( 2, static_cast<int>(log2f(indices_size-1) ) );
        int left_size = indices_size - v/2;

        if( left_size > v )
        {
            left_size = v;
        }
        int right_size = indices_size - left_size;

        float split_value = ( V.elements[current_dim+static_cast<int>(sorted_indices[current_dim].elements[left_size-1])*V.dim ] + V.elements[current_dim+static_cast<int>(sorted_indices[current_dim].elements[left_size] ) * V.dim] ) /2.0;

        kd_tree.elements[ position ] = split_value;

        LBPointArray *sorted_indices_left = (LBPointArray*)malloc( 3*sizeof(LBPointArray) );
        LBPointArray *sorted_indices_right = (LBPointArray*)malloc( 3*sizeof(LBPointArray) );

        for( int i=0; i<max_dim; i++ )
        {

            sorted_indices_left[i].width = left_size;
            sorted_indices_left[i].dim = 1;
            sorted_indices_left[i].elements = (float*)malloc( (left_size+1) *sizeof(float) );

            sorted_indices_right[i].width = right_size;
            sorted_indices_right[i].dim = 1;
            sorted_indices_right[i].elements = (float*)malloc( (right_size+1) * sizeof(float) );

            if( i == current_dim ){
                splitPointArray( sorted_indices[i], sorted_indices_left[i], sorted_indices_right[i]);
            } else {
                splitPointArrayWithValue(V, sorted_indices[i], sorted_indices_left[i], sorted_indices_right[i], current_dim, split_value);
            }

        }

        //threadeeen
        int next_dim = (current_dim+1)%max_dim;
        
        // thread pool when split
        // on depth 2 you can use 8 threads
        if(current_depth == st_depth_threads )
        {
            pool->push(generateKdTreeRecursive, V, sorted_indices_left, next_dim, max_dim, kd_tree, size, max_tree_depth, left, current_depth + 1);
            pool->push(generateKdTreeRecursive, V, sorted_indices_right, next_dim, max_dim, kd_tree, size, max_tree_depth, right, current_depth +1);
        } else {
            generateKdTreeRecursive(0, V, sorted_indices_left, (current_dim+1)%max_dim, max_dim, kd_tree, size, max_tree_depth, left, current_depth + 1);
            generateKdTreeRecursive(0, V, sorted_indices_right, (current_dim+1)%max_dim, max_dim, kd_tree, size, max_tree_depth, right, current_depth +1);
        }

    }

    for(int i=0; i<max_dim; i++) {
        free(sorted_indices[i].elements );
    }

}


void LBKdTree::sortByDim(LBPointArray& V, int dim, LBPointArray& indices, LBPointArray& values) {

    naturalMergeSort(V, dim, indices, values);

}

void LBKdTree::naturalMergeSort(LBPointArray& in, int dim, LBPointArray& indices, LBPointArray& m, int limit) {

    copyDimensionToPointArray(in, dim, m);

    int m_elements = m.width * m.dim;

    int slide_buffer_size = int(m_elements-0.5);
    int* slide_buffer = (int*) malloc(slide_buffer_size * sizeof(int));


    //create RUNS
    int num_slides = 1;
    slide_buffer[0] = 0;
    for(int i=1; i < slide_buffer_size+1; i++)
    {
        if(m.elements[i] < m.elements[i-1])
        {
            slide_buffer[num_slides] = i;
            num_slides++;
        }

    }
    slide_buffer[num_slides] = m_elements;
    slide_buffer_size = num_slides+1;


    //sort
    int count = 0;
    int current_limit = -1;
    while(num_slides > 1)
    {
        if(num_slides > 2){
            current_limit = limit;
        }

        int i;

        for(i=2;i<int(num_slides+1);i+=2)
        {

            mergeHostWithIndices(m.elements, indices.elements , slide_buffer[i-2], slide_buffer[i-1]-1, slide_buffer[i-1], slide_buffer[i]-1, current_limit);


            slide_buffer[i/2]= slide_buffer[i];
        }

        if(num_slides%2 == 1){
            slide_buffer[(num_slides+1)/2] = slide_buffer[num_slides];
        }

        count ++;
        num_slides = int(num_slides/2.0+0.5);

    }

    free(slide_buffer);
}

void LBKdTree::mergeHostWithIndices(float* a, float* b, int i1, int j1, int i2, int j2, int limit) {

    int limit_end = limit;

    float* temp = (float*) malloc((j2-i1+1) * sizeof(float));  //array used for merging
    int* temp_indices = (int*) malloc((j2-i1+1) * sizeof(int));  //array used for merging

    int i,j,k;
    i=i1;    //beginning of the first list
    j=i2;    //beginning of the second list
    k=0;

    int counter = 0;

    while( i<=j1 && j<=j2 && limit!=0 )    //while elements in both lists
    {
        counter ++;
        limit--;
        if(a[i]<a[j]){
            temp_indices[k] = b[i];
            temp[k++]=a[i++];

        }else{
            temp_indices[k] = b[j];
            temp[k++]=a[j++];
        }
    }

    while(i <= j1 && limit != 0) //copy remaining elements of the first list
    {
        temp_indices[k] = b[i];
        temp[k++]=a[i++];
    }

    while(j <= j2 && limit!=0 ) {   //copy remaining elements of the second list
        temp_indices[k] = b[j];
        temp[k++]=a[j++];
    }

    //Transfer elements from temp[] back to a[]
    for(i=i1,j=0;i<=j2 && limit_end!=0 ;i++,j++,limit_end--)
    {
        b[i] = temp_indices[j];
        if(b[i] < 0){
            printf("THERE IS SOMETHING WRONG\n");
        }
        a[i] = temp[j];
    }

    free(temp_indices);
    free(temp);
}
