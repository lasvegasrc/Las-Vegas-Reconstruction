#ifndef OVERLAPPING_KD_TREE_HPP
#define OVERLAPPING_KD_TREE_HPP

#include "OverlappingNode.hpp"
#include <iostream>
#include <boost/shared_ptr.hpp>

namespace lvr{

template<typename VertexT>
class OverlappingKdTree{
    public:
        OverlappingKdTree(unsigned int max_leaf_size=10000000);

        void insert(VertexT point);

        void calculateOverlaps(double overlap_size);

        void finishConstruction();

        size_t getNumLeafs();
    
        void printKdTreeLeafs();

    private:

        void collectLeafs(boost::shared_ptr<OverlappingNode<VertexT> > node);

        void printKdTreeLeafs(boost::shared_ptr<OverlappingNode<VertexT> > node, unsigned int depth);

        boost::shared_ptr<OverlappingNode<VertexT> > m_root;

        boost::shared_ptr<OverlappingNode<VertexT> > m_iterator;

        std::vector< OverlappingNode<VertexT>* > m_leafs;

        unsigned int m_max_leaf_size;

        unsigned int m_index_counter;

        double m_overlap_size;
};

} //namespace lvr

#include "OverlappingKdTree.tcc"

#endif //OVERLAPPING_KD_TREE_HPP