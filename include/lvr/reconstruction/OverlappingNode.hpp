#ifndef OVERLAPPING_NODE_HPP
#define OVERLAPPING_NODE_HPP

#include <vector>
#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include <iostream>
#include <lvr/geometry/BoundingBox.hpp>

namespace lvr{

template<typename VertexT>
class OverlappingNode{
    public:
        OverlappingNode(unsigned int bucket_size);

        // Node functions

        void clearNodeData();

        bool isLeaf();

        void calculateOverlap(double overlap_size);

        boost::shared_ptr<OverlappingNode<VertexT> > leftChild();
        void setLeftChild(OverlappingNode<VertexT>* node );
        void setLeftChild(boost::shared_ptr<OverlappingNode<VertexT> > node);

        boost::shared_ptr<OverlappingNode<VertexT> > rightChild();
        void setRightChild(OverlappingNode<VertexT>* node );
        void setRightChild(boost::shared_ptr<OverlappingNode<VertexT> > node);

        OverlappingNode<VertexT>* getParent();
        void setParent(OverlappingNode<VertexT>* node );

        double getSplitValue();

        unsigned int getSplitDimension();

        // Leaf functions

        void insert(VertexT point);

        unsigned int getLeafSize();

        

        std::vector< std::vector<VertexT> > getNeighboursOverlap();

        boost::shared_ptr< std::vector<VertexT> > getPoints();
        
        std::vector<VertexT>& getOverlap(unsigned int lap_id);

    private:

        OverlappingNode<VertexT>* getNeighbour();

        void split();

        unsigned int getLongestSideDim();

        // Node is Node
        boost::shared_ptr<OverlappingNode<VertexT> > m_left;
        boost::shared_ptr<OverlappingNode<VertexT> > m_right;
        OverlappingNode<VertexT>* m_parent;

        double m_split_value;

        BoundingBox<VertexT> m_bounding_box;
        

        unsigned int m_bucket_size;

        unsigned int m_split_dim;

        // Node is Leaf
        // ALL POINTS
        boost::shared_ptr< std::vector<VertexT> > m_points;

        // Overlap
        std::vector< std::vector<VertexT> > m_overlaps;
        
};

} //namespace lvr

#include "OverlappingNode.tcc"

#endif //OVERLAPPING_NODE_HPP