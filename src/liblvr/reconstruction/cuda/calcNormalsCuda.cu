/*
 * This file is part of cudaNormals.
 *
 * cudaNormals is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cudaNormals.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * calcNormalsCuda.h
 *
 * @author Alexander Mock
 * @author Matthias Greshake
 */

#include <lvr/reconstruction/cuda/calcNormalsCuda.h>

/// Define Kernels

namespace lvr {

__global__ void FlipNormalsKernel(const PointArray D_V, PointArray D_Result_Normals, float x, float y, float z);

__global__ void KNNKernel(const PointArray D_V, const PointArray D_kd_tree, PointArray D_Result_Normals, int k, int method);




__global__ void FlipNormalsKernel(const PointArray D_V,
                                  PointArray D_Result_Normals,
                                  float x, float y, float z)
{
	const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if(tid < D_V.width){
		float x_dir = x - D_V.elements[tid];
		float y_dir = y - D_V.elements[D_V.width + tid];
		float z_dir = z - D_V.elements[2 * D_V.width + tid];
		
		float scalar = ( x_dir * D_Result_Normals.elements[tid] + y_dir * D_Result_Normals.elements[D_Result_Normals.width + tid] + z_dir * D_Result_Normals.elements[2 * D_Result_Normals.width + tid] );
		
		// gegebenfalls < durch > ersetzen
		if(scalar < 0)
		{
			D_Result_Normals.elements[tid] = -D_Result_Normals.elements[tid];
			D_Result_Normals.elements[D_Result_Normals.width + tid] = -D_Result_Normals.elements[D_Result_Normals.width + tid];
			D_Result_Normals.elements[2 * D_Result_Normals.width + tid] = -D_Result_Normals.elements[2 * D_Result_Normals.width + tid];
		}
	}
}

// Get a matrix element
__device__ int GetKdTreePosition(const PointArray& D_kd_tree, float x, float y, float z)
{
	int pos = 0;
	int current_dim = 0;
	
	while(pos*2+1 < D_kd_tree.width)
	{
		
		if(current_dim == 0)
		{
			if(x <= D_kd_tree.elements[pos] )
			{
				pos = pos*2+1;
			} else {
				pos = pos*2+2;
			}
			
			current_dim += 1;
			
		} else if(current_dim == 1) {
			
			if(y <= D_kd_tree.elements[pos] ){
				pos = pos*2+1;
			}else{
				pos = pos*2+2;
			}
			
			current_dim +=1;
		} else {
			if(z <= D_kd_tree.elements[pos] ){
				pos = pos*2+1;
			}else{
				pos = pos*2+2;
			}
			
			current_dim = 0;
		}
		
	}
	
    return pos;
}

__device__ float SearchQueryPoint(const PointArray& D_kd_tree, float x, float y, float z)
{
	return D_kd_tree.elements[GetKdTreePosition(D_kd_tree, x, y, z)];
}


__device__ void calculateNormalRansa2(float* nn_vecs,
                                      int k,
                                      int max_iterations,
                                      float& x, float& y, float& z)
{
	float min_dist = FLT_MAX;
	int iterations = 0;
	
	for(int i=3; i<k*3; i+=3){
		//~ printf("%f %f %f\n", last_vec[0], last_vec[1], last_vec[2]);
		
		int j = (i + int(k/3) * 3) % (k * 3);
		
		float n_x = nn_vecs[j+1]*nn_vecs[i+2] - nn_vecs[j+2]*nn_vecs[i+1];
		float n_y = nn_vecs[j+2]*nn_vecs[i+0] - nn_vecs[j+0]*nn_vecs[i+2];
		float n_z = nn_vecs[j+0]*nn_vecs[i+1] - nn_vecs[j+1]*nn_vecs[i+0];
		
		float norm = sqrtf( n_x*n_x + n_y*n_y + n_z*n_z );
		
		
		if( norm != 0.0 ){
			
			float norm_inv = 1.0/norm;
			
			n_x = n_x * norm_inv;
			n_y = n_y * norm_inv;
			n_z = n_z * norm_inv;
			
			float cum_dist = 0.0;
			for(int j=0; j<k*3; j+=3){
				cum_dist += abs(n_x * nn_vecs[j] + n_y * nn_vecs[j+1] + n_z * nn_vecs[j+2]);
			}
			
			if(cum_dist < min_dist) {
				
				iterations = 0;
				min_dist = cum_dist;
				x = n_x;
				y = n_y;
				z = n_z;
				
			} else if(iterations < max_iterations) {
				
				iterations ++;
			
			}else{
				
				return;
				
			}
		}
		
	}
}

__device__ void calculateNormalRansa(float* nn_vecs,
                                     int k,
                                     int max_iterations,
                                     float& x, float& y, float& z)
{
	
	float * last_vec = (float*)malloc(3 * sizeof(float) );
	last_vec[0] = nn_vecs[0];
	last_vec[1] = nn_vecs[1];
	last_vec[2] = nn_vecs[2];
	float min_dist = FLT_MAX;
	int iterations = 0;
	// nearest neighbors in nn!!
	// what now? 
	// PCA?
	// minimize plane error:
	
	
	for(int i=3; i<k*3; i+=3){
		// cross product

		//~ printf("%f %f %f\n", last_vec[0], last_vec[1], last_vec[2]);
		float n_x = last_vec[1]*nn_vecs[i+2] - last_vec[2]*nn_vecs[i+1];
		float n_y = last_vec[2]*nn_vecs[i+0] - last_vec[0]*nn_vecs[i+2];
		float n_z = last_vec[0]*nn_vecs[i+1] - last_vec[1]*nn_vecs[i+0];
		

		float norm = sqrtf( n_x*n_x + n_y*n_y + n_z*n_z );
		
		if( norm == 0.0){
				
			last_vec[0] = nn_vecs[i+0];
			last_vec[1] = nn_vecs[i+1];
			last_vec[2] = nn_vecs[i+2];
			continue;
			
                }

		float norm_inv = 1.0/norm;
		//~ float norm = n_x*n_x + n_y*n_y + n_z*n_z ;
		n_x = n_x * norm_inv;
		n_y = n_y * norm_inv;
		n_z = n_z * norm_inv;
		//~ printf("%f %f %f\n",n_x,n_y,n_z);
		
		float cum_dist = 0.0;
		for(int j=0; j<k*3; j+=3){
			cum_dist += abs(n_x * nn_vecs[j] + n_y * nn_vecs[j+1] + n_z * nn_vecs[j+2]);
		}
		
		if(cum_dist < min_dist){
			iterations = 0;
			min_dist = cum_dist;
			x = n_x;
			y = n_y;
			z = n_z;
			//~ printf("%f %f %f\n",x,y,z);
		}else{
			iterations +=1;
		}
		
		last_vec[0] = nn_vecs[i+0];
		last_vec[1] = nn_vecs[i+1];
		last_vec[2] = nn_vecs[i+2];
		
		if(iterations > max_iterations){
			break;
		}
	}
	
	//instead of minimize plane error:
	// take normal with maximum of inliers (RANSAC like)
	
	free(last_vec);
}

__device__ void calculateNormalPCA(float* nn_vecs, int k, float& n_x, float& n_y, float& n_z)
{
	
	// ilikebigbits.com/blog/2015/3/2/plane-from-points
	
	
	//x
	float xx = 0.0;
	float xy = 0.0;
	float xz = 0.0;
	
	//y
	float yy = 0.0;
	float yz = 0.0;
	
	//z
	float zz = 0.0;
	
	for(int i=0; i<k; i++)
	{
		float rx = nn_vecs[i*3+0];
		float ry = nn_vecs[i*3+1];
		float rz = nn_vecs[i*3+2];
		
		xx += rx * rx;
		xy += rx * ry;
		xz += rx * rz;
		yy += ry * ry;
		yz += ry * rz;
		zz += rz * rz;
	}
	
	//determinante? 
	float det_x = yy * zz - yz * yz;
	float det_y = xx * zz - xz * xz;
	float det_z = xx * yy - xy * xy;
	
	float dir_x;
	float dir_y;
	float dir_z;
	// det X biggest
	if( det_x >= det_y && det_x >= det_z){
		
		if(det_x <= 0.0){
			//not a plane
		}
		
		dir_x = 1.0;
		dir_y = (xz * yz - xy * zz) / det_x;
		dir_z = (xy * yz - xz * yy) / det_x;
	} //det Y biggest
	else if( det_y >= det_x && det_y >= det_z){
		
		if(det_y <= 0.0){
			// not a plane
		}
		
		dir_x = (yz * xz - xy * zz) / det_y;
		dir_y = 1.0;
		dir_z = (xy * xz - yz * xx) / det_y;
	} // det Z biggest
	else{
		if(det_z <= 0.0){
			// not a plane
		}
		
		dir_x = (yz * xy - xz * yy ) / det_z;
		dir_y = (xz * xy - yz * xx ) / det_z;
		dir_z = 1.0;
	}
	
	float invnorm = 1/sqrtf( dir_x * dir_x + dir_y * dir_y + dir_z * dir_z );
	
	n_x = dir_x * invnorm;
	n_y = dir_y * invnorm;
	n_z = dir_z * invnorm;
	
}

__device__ void switchNeighbor(float* nn_vecs, int k, float v_x, float v_y, float v_z)
{
	
	if( ( v_x==0.0 || v_x==-0.0 ) &&
			( v_y==0.0 || v_y==-0.0 ) && 
			( v_z==0.0 || v_z==-0.0 ) )
	{
		return;
	}
	
	for(int i=0; i<k*3; i+=3){
		if( ( nn_vecs[i]==0.0 || nn_vecs[i]==-0.0 ) &&
			( nn_vecs[i+1]==0.0 || nn_vecs[i+1]==-0.0 ) && 
			( nn_vecs[i+2]==0.0 || nn_vecs[i+2]==-0.0 ) )
		{
			nn_vecs[i] = v_x;
			nn_vecs[i+1] = v_y;
			nn_vecs[i+2] = v_z;
		} else {
			float dist_old = nn_vecs[i]*nn_vecs[i] + nn_vecs[i+1]*nn_vecs[i+1] + nn_vecs[i+2]*nn_vecs[i+2];
			float dist_new = v_x*v_x + v_y*v_y + v_z*v_z;
			if(dist_new < dist_old){
				nn_vecs[i] = v_x;
				nn_vecs[i+1] = v_y;
				nn_vecs[i+2] = v_z;
			}
		}
	}
}

__device__ void getNearestNeighbors(const PointArray& D_V,
                                    const PointArray& D_kd_tree,
                                    int k,
                                    int subtree_pos,
                                    int pos,
                                    int pos_value,
                                    float* nn_vecs)
{
	
    int iterator = subtree_pos;
    int max_nodes = 1;
    bool leaf_reached = false;
    int i_nn = 0;

    int query_index = pos_value * D_V.dim;

    float query_x = D_V.elements[ query_index ];
    float query_y = D_V.elements[ query_index + 1 ];
    float query_z = D_V.elements[ query_index + 2 ];

    // like width search
    // go kd-tree up until max_nodes(leaf_nodes of subtree) bigger than needed nodes k
    // iterator = iterator * 2 + 1 -> go to
    for( ; iterator < D_kd_tree.width; iterator = iterator * 2 + 1, max_nodes *= 2)
    {
    // collect nodes from current height
        for( int i=0; i < max_nodes && iterator + i < D_kd_tree.width; i++)
        {
            int current_pos = iterator + i;
            int leaf_value  = (int)(D_kd_tree.elements[ current_pos ] + 0.5 );

            if( leaf_reached && i_nn <= k*3 )
            {
                if( leaf_value != pos_value && leaf_value < D_V.width )
                {
                    int curr_nn_index = leaf_value * D_V.dim;

                    float nn_x = D_V.elements[ curr_nn_index ] - query_x;
                    float nn_y = D_V.elements[ curr_nn_index + 1 ] - query_y;
                    float nn_z = D_V.elements[ curr_nn_index + 2 ] - query_z;

                    if(nn_x != 0.0 || nn_y != 0.0 || nn_z != 0.0)
                    {
                        nn_vecs[ i_nn ]     = nn_x;
                        nn_vecs[ i_nn + 1 ] = nn_y;
                        nn_vecs[ i_nn + 2 ] = nn_z;
                        i_nn += 3;
                    }
                }
            } else if( current_pos * 2 + 1 >= D_kd_tree.width ) {

                int curr_nn_index = leaf_value * D_V.dim;
                //first leaf reached
                leaf_reached = true;
                if( leaf_value != pos_value && i_nn <= k*3 )
                {
                    nn_vecs[i_nn]   = D_V.elements[ curr_nn_index ] - query_x;
                    nn_vecs[i_nn+1] = D_V.elements[ curr_nn_index + 1 ] - query_y;
                    nn_vecs[i_nn+2] = D_V.elements[ curr_nn_index + 2] - query_z;
                    i_nn += 3;
                }
            }
        }
    }
}

__device__ bool checkLinearNeighborHood(const PointArray& D_V, const PointArray& D_kd_tree, int pos, int k){
	
	int number_true = 0;
	int * split_positions = (int*)malloc(6*sizeof(int));
	split_positions[0] = (int)((pos  - 1) / 2);
	
	for(int i=1; i<6; i++){
		split_positions[i] = (int)((split_positions[i-1]  - 1) / 2);
		
        }
	
	// check x
	for(int i=0;i<3;i++)
	{
		
		if(split_positions[i+3] > 0 )
                {
			if(D_kd_tree.elements[split_positions[i+3] ] != D_kd_tree.elements[split_positions[i] ] )
			{	
				number_true += 1;
			}
		}else{
			number_true += 1;
		}
	}
	
	
	free(split_positions);
	
	if(number_true >= 2){
		return false;
	}else{
		return true;
	}
}

__device__ void calculateNormalFromSubtree(const PointArray& D_V, const PointArray& D_kd_tree, int pos, int k, float& x, float& y, float& z, int method )
{
        //~
        //~  Step 1: get upper node
        //~  Step 2: get child nodes != query node
        //~  Step 3: calculate normals
        //~

        int pos_value = (int)(D_kd_tree.elements[pos]+0.5);

        int subtree_pos = pos;
        int i;
        for(i=1; i<(k+1) && subtree_pos>0; i*=2) {
                subtree_pos = (int)((subtree_pos  - 1) / 2);
        }
        //~ printf("subtree_pos: %d\n",subtree_pos);

        // k+1 FIX
        float * nn_vecs = (float*)malloc(3*(k+1)*sizeof(float));


        getNearestNeighbors(D_V, D_kd_tree, k, subtree_pos, pos, pos_value, nn_vecs);

        if(method == 0){
                //PCA
                calculateNormalPCA(nn_vecs, k, x, y, z);
        } else if(method == 1) {
                //RANSAC
                calculateNormalRansa2(nn_vecs, k, 8, x, y, z);
        }

        free(nn_vecs);
	
} 

//distance function without transformation
__global__ void KNNKernel(const PointArray D_V, const PointArray D_kd_tree, PointArray D_Result_Normals, int k, int method)
{
	const unsigned int tid = blockDim.x * blockIdx.x + threadIdx.x;
	
	if(tid < D_V.width){
		
		int pos = GetKdTreePosition(D_kd_tree, D_V.elements[tid * D_V.dim], D_V.elements[tid * D_V.dim + 1], D_V.elements[tid * D_V.dim +2] );
		
		float result_x = D_Result_Normals.elements[tid * D_Result_Normals.dim ];
		float result_y = D_Result_Normals.elements[tid * D_Result_Normals.dim + 1 ];
		float result_z = D_Result_Normals.elements[tid * D_Result_Normals.dim + 2 ];
		
                calculateNormalFromSubtree(D_V, D_kd_tree, pos, k, result_x, result_y, result_z, method);
		
		D_Result_Normals.elements[tid * D_Result_Normals.dim ] = result_x;
		D_Result_Normals.elements[tid * D_Result_Normals.dim + 1 ] = result_y;
		D_Result_Normals.elements[tid * D_Result_Normals.dim + 2 ] = result_z;
		
	}
	
}

/// HOST FUNCTIONS ///

void CalcNormalsCuda::init(){
	// set default k
	this->m_k = 50;
	
	// set default flippoint
	this->m_vx = 1000000.0;
	this->m_vy = 1000000.0;
	this->m_vz = 1000000.0;
	
	this->m_calc_method = 0;

    this->m_b_factor = 16;
}

CalcNormalsCuda::CalcNormalsCuda(PointArray& points)
{
	this->init();
	
	this->getCudaInformation();
	
	this->V.dim = points.dim;
	
	this->V.width = points.width;
	
	mallocPointArray(V);
	
	for(int i = 0; i<points.width*points.dim; i++)
	{
		
		this->V.elements[i] = points.elements[i];
		
	}
	
	this->initKdTree();
	
}

CalcNormalsCuda::CalcNormalsCuda(floatArr& points, size_t num_points, size_t dim){
	
	this->init();
	
	this->getCudaInformation();
	
	this->V.dim = static_cast<int>(dim);
	
	this->V.width = static_cast<int>(num_points);
	
	mallocPointArray(V);
	
	this->V.elements = points.get();
	
	
	this->initKdTree();
	
}

void CalcNormalsCuda::getCudaInformation()
{
	
	m_mps = 0;
	m_threads_per_mp = 0;
	m_threads_per_block = 0;
	m_size_thread_block = new int(3);
	m_size_grid = new int(3);
	m_device_global_memory = 0;
	
	
	cudaSetDevice(0);
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, 0);
    
    
        m_mps = deviceProp.multiProcessorCount;
        m_threads_per_mp = deviceProp.maxThreadsPerMultiProcessor;
        m_threads_per_block = deviceProp.maxThreadsPerBlock;
        m_size_thread_block[0] = deviceProp.maxThreadsDim[0];
        m_size_thread_block[1] = deviceProp.maxThreadsDim[1];
        m_size_thread_block[2] = deviceProp.maxThreadsDim[2];
        m_size_grid[0] = deviceProp.maxGridSize[0];
        m_size_grid[1] = deviceProp.maxGridSize[1];
        m_size_grid[2] = deviceProp.maxGridSize[2];
        m_device_global_memory = (unsigned long long) deviceProp.totalGlobalMem;
    
}

void CalcNormalsCuda::getNormals(PointArray& output_normals)
{
	
	output_normals.dim = this->Result_Normals.dim;
	output_normals.width = this->Result_Normals.width;
	//output_normals.elements = (float*)malloc( this->Result_Normals.dim * this->Result_Normals.width * sizeof(float) ) ;
	
	for(int i = 0; i< this->Result_Normals.dim * this->Result_Normals.width; i++)
	{	
		output_normals.elements[i] = this->Result_Normals.elements[i];
	}
	
}

void CalcNormalsCuda::getNormals(floatArr output_normals)
{
	for(int i = 0; i< this->Result_Normals.dim * this->Result_Normals.width; i++)
	{	
		output_normals[i] = this->Result_Normals.elements[i];
	}
}

void CalcNormalsCuda::mallocPointArray(PointArray& m) {

	m.elements = (float*)malloc(m.width * m.dim * sizeof(float));

}

void CalcNormalsCuda::generateHostPointArray(PointArray& m, int width, int dim)
{
	
	m.dim = dim;
	m.width = width;
	m.elements = (float*)malloc(m.width * m.dim * sizeof(float) );
	
}

void CalcNormalsCuda::generateDevicePointArray(PointArray& D_m, int width, int dim) {
	
        D_m.width = width;
        D_m.dim = dim;
        size_t size = D_m.width * D_m.dim * sizeof(float);
        HANDLE_ERROR( cudaMalloc(&D_m.elements, size) );
    
}

void CalcNormalsCuda::copyToDevicePointArray(PointArray& m, PointArray& D_m) {
	
	size_t size = m.width * m.dim * sizeof(float);
        HANDLE_ERROR( cudaMemcpy(D_m.elements, m.elements, size, cudaMemcpyHostToDevice) );

}

void CalcNormalsCuda::copyToHostPointArray(PointArray& D_m, PointArray& m) {
	
	size_t size = m.width * m.dim * sizeof(float);
        HANDLE_ERROR( cudaMemcpy(m.elements, D_m.elements, size, cudaMemcpyDeviceToHost) );
}

void CalcNormalsCuda::fillPointArrayWithSequence(PointArray& m) {

        for(int i=0;i<m.width*m.dim;i++)
        {
                m.elements[i] = i;
        }

}  

void CalcNormalsCuda::copyDimensionToPointArray(PointArray& in, int dim, PointArray& out) {

        for(int i = 0; i<out.width; i++)
        {
		out.elements[i] = in.elements[i * in.dim + dim];
	}
}

void CalcNormalsCuda::copyVectorInterval(PointArray& in, int start, int end, PointArray& out) {

        for(int i=0; i < (end-start); i++)
        {
		out.elements[i] = in.elements[i + start];
	}
}

void CalcNormalsCuda::mergeHostWithIndices(float* a, float* b, int i1, int j1, int i2, int j2, int limit) {

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
                if( a[i] < a[j] ) {
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


void CalcNormalsCuda::naturalMergeSort(PointArray& in, int dim, PointArray& indices, PointArray& m, int limit) {
	
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

            mergeHostWithIndices( m.elements, indices.elements , slide_buffer[i-2], slide_buffer[i-1]-1, slide_buffer[i-1], slide_buffer[i]-1, current_limit );
			
			
            slide_buffer[i/2] = slide_buffer[i];
		}
		
		if(num_slides%2 == 1){
            slide_buffer[ (num_slides+1)/2 ] = slide_buffer[num_slides];
		}
		
		count ++;
		num_slides = int(num_slides/2.0+0.5);
		
	}
	
	free(slide_buffer);
}

void CalcNormalsCuda::sortByDim(PointArray& V, int dim, PointArray& indices, PointArray& values) {

        naturalMergeSort(V, dim, indices, values);

}

void CalcNormalsCuda::splitPointArray(PointArray& I, PointArray& I_L, PointArray& I_R) {
	
	int i=0;
	for(; i < I_L.width * I_L.dim; i++){
		I_L.elements[i] = I.elements[i];
	}
	int j=0;
	for(; i<I.width*I.dim && j<I_R.width*I_R.dim; i++, j++){
		I_R.elements[j] = I.elements[i];
	}
	
}

void CalcNormalsCuda::splitPointArrayWithValue(PointArray& V, PointArray& I, PointArray& I_L, PointArray& I_R, int current_dim, float value) {

    int i_l = 0;
	int i_r = 0;
	
	for(int i=0; i<I.width; i++)
	{
		float current_value = V.elements[static_cast<int>(I.elements[i] + 0.5) * V.dim + current_dim ];
		//~ printf("curr val: %f\n", current_value);
		if(current_value <= value && I_L.width > i_l ){
			//~ printf("add to left: %f with value %f\n", I.elements[i], current_value);
			I_L.elements[i_l++] = I.elements[i];
		}else if(current_value >= value && I_R.width > i_r){
			//~ printf("add to right: %f with value %f\n", I.elements[i], current_value);
			I_R.elements[i_r++] = I.elements[i];
		}else {
			if(i_r<I_R.width){
				I_R.elements[i_r++] = I.elements[i];
			}else if(i_l<I_L.width){
				I_L.elements[i_l++] = I.elements[i];
			}
		}
	}
		
}


void CalcNormalsCuda::generateKdTreeRecursive(PointArray& V, PointArray* sorted_indices, int current_dim, int max_dim, PointArray& kd_tree, int size, int max_tree_depth, int position) {
	
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
		
		struct PointArray sorted_indices_left[max_dim];
		struct PointArray sorted_indices_right[max_dim];
		
		// alloc new memory
        for( int i=0; i<max_dim; i++ )
        {
			//memory corruption when malloc 
			
			sorted_indices_left[i].width = left_size;
			sorted_indices_left[i].dim = 1;
			sorted_indices_left[i].elements = (float*)malloc( (left_size+1) *sizeof(float) );
			
			sorted_indices_right[i].width = right_size;
			sorted_indices_right[i].dim = 1;
			sorted_indices_right[i].elements = (float*)malloc( (right_size+1) * sizeof(float) );
			
            if( i == current_dim ) {
				splitPointArray( sorted_indices[i], sorted_indices_left[i], sorted_indices_right[i]);
			}else{
				splitPointArrayWithValue(V, sorted_indices[i], sorted_indices_left[i], sorted_indices_right[i], current_dim, split_value);
			}
			
		}
		
		generateKdTreeRecursive(V, sorted_indices_left, (current_dim+1)%max_dim, max_dim, kd_tree, size, max_tree_depth, left);
		generateKdTreeRecursive(V, sorted_indices_right, (current_dim+1)%max_dim, max_dim, kd_tree, size, max_tree_depth, right);		
	
		
		// alloc new memory
        for(int i=0; i<max_dim; i++)
        {
            free( sorted_indices_left[i].elements );
            free( sorted_indices_right[i].elements );
		}
	}
	
	
}

void CalcNormalsCuda::generateKdTreeArray(PointArray& V, PointArray* sorted_indices, int max_dim, PointArray& kd_tree) {

    int size;
	int max_tree_depth;
	
	max_tree_depth = static_cast<int>( log2f(V.width - 1 ) + 2.0 ) ;

	if (V.width == 1)
	{
		max_tree_depth = 1;
	}
	
	size = V.width * 2 - 1;
	
        generateHostPointArray(kd_tree, size, 1);
	
	//start real generate
	generateKdTreeRecursive(V, sorted_indices, 0, max_dim, kd_tree, size, max_tree_depth, 0);
	
}

void CalcNormalsCuda::GPU_NN(PointArray& D_V, PointArray& D_kd_tree, PointArray& D_Result_Normals, PointArray& Result_Normals) {
	
    int threadsPerBlock = this->m_threads_per_block;
	int blocksPerGrid = (D_V.width + threadsPerBlock-1)/threadsPerBlock;

    // kNN-search and Normal calculation
    KNNKernel<<<blocksPerGrid * this->m_b_factor, threadsPerBlock / this->m_b_factor >>>(D_V, D_kd_tree, D_Result_Normals, this->m_k, this->m_calc_method);
    cudaDeviceSynchronize();

    // Flip normals to view point
    FlipNormalsKernel<<<blocksPerGrid * this->m_b_factor, threadsPerBlock / this->m_b_factor >>>(D_V, D_Result_Normals, this->m_vx, this->m_vy, this->m_vz);
	cudaDeviceSynchronize();

    //TODO: Interpolate
	
    size_t size = Result_Normals.width * Result_Normals.dim * sizeof(float);

    HANDLE_ERROR( cudaMemcpy(Result_Normals.elements, D_Result_Normals.elements, size, cudaMemcpyDeviceToHost ) );

}

void CalcNormalsCuda::initKdTree() {
	
	//~ struct Matrix test;
	struct PointArray indices_sorted[this->V.dim];
	struct PointArray values_sorted[this->V.dim];
	
	for(int i=0; i < this->V.dim; i++)
	{
        generateHostPointArray(indices_sorted[i], V.width, 1);
		
        generateHostPointArray(values_sorted[i], V.width, 1);
		
		fillPointArrayWithSequence(indices_sorted[i]);
		
		sortByDim( this->V, i, indices_sorted[i] , values_sorted[i]);
	}

	generateKdTreeArray(V, indices_sorted, this->V.dim, this->kd_tree);
	
	for(int i=0; i<V.dim;i++)
	{
		free(indices_sorted[i].elements);
		free(values_sorted[i].elements);
	}

}

void CalcNormalsCuda::setK(int k) {

    this->m_k = k;

}

void CalcNormalsCuda::setFlippoint(float v_x, float v_y, float v_z) {

    this->m_vx = v_x;
	this->m_vy = v_y;
	this->m_vz = v_z;

}

void CalcNormalsCuda::setMethod(std::string method) {

    if( strcmp( method.c_str(), "PCA") == 0 ){
		this->m_calc_method = 0;
	} else if( strcmp( method.c_str(), "RANSAC") == 0){
		this->m_calc_method = 1;
	} else {
		printf("WARNING: Normal Calculation Method is not implemented\n");
	}

}

void CalcNormalsCuda::setBlockSizeFactor(int b_factor) {
    this->m_b_factor = b_factor;
}

void CalcNormalsCuda::printSettings() {

	printf("	Nearest Neighbors = %d\n",this->m_k);
	
	printf("	Flip point = (%f, %f, %f)\n", this->m_vx, this->m_vy, this->m_vz);
	
	switch(this->m_calc_method){
		case 0:
			printf("	Method = 'PCA'\n");
			break;
		case 1:
			printf("	Method = 'RANSAC'\n");
			break;
	}
	
	printf("\n");

}

void CalcNormalsCuda::start() {
	
	generateHostPointArray( this->Result_Normals, V.width, V.dim);
	
	PointArray D_V, D_kd_tree, D_Result_Normals;
	generateDevicePointArray( D_V, this->V.width, this->V.dim );
	generateDevicePointArray( D_kd_tree, this->kd_tree.width, this->kd_tree.dim);
	generateDevicePointArray( D_Result_Normals, this->Result_Normals.width, this->Result_Normals.dim);
	
	//COPY STUFF
	copyToDevicePointArray( V, D_V);
	copyToDevicePointArray( this->kd_tree, D_kd_tree);
	
	//Cuda Kernels
	GPU_NN(D_V, D_kd_tree, D_Result_Normals, this->Result_Normals);
	
	cudaFree(D_V.elements);
	cudaFree(D_kd_tree.elements);
	cudaFree(D_Result_Normals.elements);
		
}

CalcNormalsCuda::~CalcNormalsCuda() {

    // clearn up resulting normals and kd_tree
    // Pointcloud has to be cleaned up by user

	if(this->Result_Normals.width > 0){
		free(Result_Normals.elements);
	}
	
	if(this->kd_tree.width > 0){
		free(this->kd_tree.elements);
	}
}

} /* namespace lvr */
