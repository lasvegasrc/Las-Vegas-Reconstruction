#ifndef __LBKDTREE_HPP
#define __LBKDTREE_HPP

#include "LBPointArray.hpp"
#include "ctpl.h"

#include <stdlib.h>
#include <math.h>



/**
 * @brief The LBKdTree class implements a left-balanced array-based index kd-tree.
 *          Left-Balanced: minimum memory
 *          Array-Based: Good for GPU - Usage
 */
class LBKdTree {
public:

    LBKdTree( LBPointArray& vertices , int num_threads=8);

    void generateKdTree( LBPointArray& vertices );

    LBPointArray getKdTreeArray();

    

private:

    
    void generateKdTreeArray(LBPointArray& V, LBPointArray* sorted_indices, int max_dim, LBPointArray& kd_tree);

    //void generateAndSort(int id, LBPointArray& vertices, LBPointArray* indices_sorted, LBPointArray* values_sorted, int dim);

    void sortByDim(LBPointArray& V, int dim, LBPointArray& indices, LBPointArray& values);

    void naturalMergeSort(LBPointArray& in, int dim, LBPointArray& indices,  LBPointArray& m, int limit=-1);

    void mergeHostWithIndices(float* a, float* b, int i1, int j1, int i2, int j2, int limit=-1);

    LBPointArray kd_tree;
    

    // Static member

    static int st_num_threads;
    static int st_depth_threads;

    static ctpl::thread_pool *pool;

    static void generateKdTreeRecursive(int id, LBPointArray& V, LBPointArray* sorted_indices, int current_dim, int max_dim, LBPointArray& kd_tree, int size, int max_tree_depth, int position, int current_depth);

    static void test(int id, LBPointArray* sorted_indices);
    

};


#endif // !__LBKDTREE_HPP
