#ifndef __POINTARRAY_HPP
#define __POINTARRAY_HPP

#include <stdlib.h>
#include <iostream>



struct LBPointArray {
    int width;
    int dim;
    float* elements;
};

// static helper methods

static void mallocPointArray(LBPointArray& m) {

    m.elements = (float*)malloc(m.width * m.dim * sizeof(float));

}

static void generatePointArray(LBPointArray& m, int width, int dim)
{

    m.dim = dim;
    m.width = width;
    m.elements = (float*)malloc(m.width * m.dim * sizeof(float) );

}

static void generatePointArray(int id, LBPointArray& m, int width, int dim)
{

    m.dim = dim;
    m.width = width;
    m.elements = (float*)malloc(m.width * m.dim * sizeof(float) );

}

static void fillPointArrayWithSequence(LBPointArray& m) {

    for(int i=0;i<m.width*m.dim;i++)
    {
        m.elements[i] = i;
    }

}

// Pool function
static void fillPointArrayWithSequence(int id, LBPointArray& m) {

    for(int i=0;i<m.width*m.dim;i++)
    {
        m.elements[i] = i;
    }

}

static void copyVectorInterval(LBPointArray& in, int start, int end, LBPointArray& out) {

    for(int i=0; i < (end-start); i++)
    {
        out.elements[i] = in.elements[i + start];
    }
}

static void copyDimensionToPointArray(LBPointArray& in, int dim, LBPointArray& out) {

    for(int i = 0; i<out.width; i++)
    {
        out.elements[i] = in.elements[i * in.dim + dim];
    }
}

static void splitPointArray(LBPointArray& I, LBPointArray& I_L, LBPointArray& I_R) {
    
    int i=0;
    for(; i < I_L.width * I_L.dim; i++){
        I_L.elements[i] = I.elements[i];
    }
    int j=0;
    for(; i<I.width*I.dim && j<I_R.width*I_R.dim; i++, j++){
        I_R.elements[j] = I.elements[i];
    }

}

static void splitPointArrayWithValue(LBPointArray& V, LBPointArray& I, LBPointArray& I_L, LBPointArray& I_R, int current_dim, float value) {

    int i_l = 0;
    int i_r = 0;

    for(int i=0; i<I.width; i++)
    {
        float current_value = V.elements[static_cast<int>(I.elements[i] + 0.5) * V.dim + current_dim ];
        //~ printf("curr val: %f\n", current_value);
        if(current_value <= value && I_L.width > i_l ){
            //~ printf("add to left: %f with value %f\n", I.elements[i], current_value);
            I_L.elements[i_l++] = I.elements[i];
        } else if(current_value >= value && I_R.width > i_r){
            //~ printf("add to right: %f with value %f\n", I.elements[i], current_value);
            I_R.elements[i_r++] = I.elements[i];
        } else {
            if(i_r<I_R.width){
                I_R.elements[i_r++] = I.elements[i];
            }else if(i_l<I_L.width){
                I_L.elements[i_l++] = I.elements[i];
            }
        }
    }

}

// SORT FUNCTIONS THREADED


static void mergeHostWithIndices(float* a, float* b, int i1, int j1, int i2, int j2, int limit) {

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

static void naturalMergeSort(LBPointArray& in, int dim, LBPointArray& indices, LBPointArray& m, int limit=-1) {

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


static void sortByDim(LBPointArray& V, int dim, LBPointArray& indices, LBPointArray& values) {

    naturalMergeSort(V, dim, indices, values);

}

static void generateAndSort(int id, LBPointArray& vertices, LBPointArray* indices_sorted, LBPointArray* values_sorted, int dim)
{
    generatePointArray( *indices_sorted, vertices.width, 1);
    generatePointArray( *values_sorted, vertices.width, 1);
    
    fillPointArrayWithSequence(*indices_sorted);

    sortByDim( vertices, dim, *indices_sorted , *values_sorted);

}

#endif // !__POINTARRAY_HPP
