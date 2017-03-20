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

#ifndef __CALCNORMALSCUDA_H
#define __CALCNORMALSCUDA_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <float.h>
#include <stdint.h>

#include <cuda_runtime.h>
#include <driver_types.h>

#include <boost/shared_array.hpp>



static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

namespace lvr {

struct PointArray {
    int width;
    int dim;
    float* elements;
};

//~ Device Functions
//~ __global__ void KNNKernel(const PointArray D_V, const PointArray D_kd_tree, PointArray D_Result_Normals, int k=50);

typedef boost::shared_array<float> floatArr;

class CalcNormalsCuda {

public:
	/**
	 * @brief Constructor.
	 *
     * @param points Input Pointcloud for kd-tree construction
	 */
	CalcNormalsCuda(PointArray& points);
	
	CalcNormalsCuda(floatArr& points, size_t num_points, size_t dim = 3);
	
	~CalcNormalsCuda();
	
	/**
	 * @brief Starts calculation the normals on GPU
	 *
	 */
	void start();
	
	/**
	 * @brief Get the resulting normals of the normal calculation. After calling "start".
	 * 
	 * @param output_normals 	PointArray as return value
	 */
	void getNormals(PointArray& output_normals);
	
	void getNormals(floatArr output_normals);
	
	/**
	 * @brief Set the number of k nearest neighbors
	 *        k-neighborhood
	 *
	 * @param k             The size of the used k-neighborhood
	 *
	 */
	void setK(int k);
	
	/**
	 * @brief Set the viewpoint to orientate the normals
	 *
	 * @param v_x 	Coordinate X axis
	 * @param v_y 	Coordinate Y axis
	 * @param v_z 	Coordinate Z axis
	 *
	 */
	void setFlippoint(float v_x, float v_y, float v_z);
	
	/**
	 * @brief Set Method for normal calculation
	 *
	 * @param method   "PCA","RANSAC"
	 *
	 */
	void setMethod(std::string method);
	
	
private:
	//~ Hostfunctions
	void init();
	
	void printSettings();
	
	void getCudaInformation();
	
	void mallocPointArray(PointArray& m);

	void generateHostPointArray(PointArray& m, int width, int dim);

	void generateDevicePointArray(PointArray& D_m, int width, int dim);

	void copyToDevicePointArray(PointArray& m, PointArray& D_m);

	void copyToHostPointArray(PointArray& D_m, PointArray& m);
	
	void fillPointArrayWithSequence(PointArray& m);

	void copyDimensionToPointArray(PointArray& in, int dim, PointArray& out);
	
	void copyVectorInterval(PointArray& in, int start, int end, PointArray& out);
	
	void naturalMergeSort(PointArray& in, int dim, PointArray& indices,  PointArray& m, int limit=-1);
	
	void mergeHostWithIndices(float* a, float* b, int i1, int j1, int i2, int j2, int limit=-1);
	
	void sortByDim(PointArray& V, int dim, PointArray& indices, PointArray& values);
	
	void splitPointArray(PointArray& I, PointArray& I_L, PointArray& I_R);
	
	void splitPointArrayWithValue(PointArray& V, PointArray& I, PointArray& I_L, PointArray& I_R, int current_dim, float value);
	
	void generateKdTreeRecursive(PointArray& V, PointArray* sorted_indices, int current_dim, int max_dim, PointArray& kd_tree, int size, int max_tree_depth, int position);
	
	void generateKdTreeArray(PointArray& V, PointArray* sorted_indices, int max_dim, PointArray& kd_tree);
	
	// Divice Function
	void GPU_NN(PointArray& D_V, PointArray& D_kd_tree, PointArray& D_Result_Normals, PointArray& Result_Normals);
	
	void initKdTree();
	
	
	
	// V->points and normals
	PointArray V, kd_tree, Result_Normals;
	float m_vx, m_vy, m_vz;
	int m_k;
	
	
	int m_calc_method;

	// Device Information
	int m_mps;
	int m_threads_per_mp;
	int m_threads_per_block;
	int* m_size_thread_block;
	int* m_size_grid;
	unsigned long long m_device_global_memory;
		
};

}

#endif // !__CALCNORMALSCUDA_H
