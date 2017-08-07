namespace lvr{

template<typename VertexT>
OverlappingKdTree<VertexT>::OverlappingKdTree(unsigned int max_leaf_size)
:m_max_leaf_size(max_leaf_size), 
m_index_counter(0),
m_root(boost::shared_ptr<OverlappingNode<VertexT> >(new OverlappingNode<VertexT>(max_leaf_size) ) ),
m_overlap_size(0.15)
{
    
}

template<typename VertexT>
void OverlappingKdTree<VertexT>::finishConstruction()
{
    boost::shared_ptr<OverlappingNode<VertexT> > tree_iterator(m_root);
    
    this->collectLeafs(tree_iterator);
}

template<typename VertexT>
void OverlappingKdTree<VertexT>::insert(VertexT point)
{
    boost::shared_ptr<OverlappingNode<VertexT> > tree_iterator(m_root);
  
    tree_iterator->insert(point);
}

template<typename VertexT>
void OverlappingKdTree<VertexT>::calculateOverlaps(double overlap_size)
{
    boost::shared_ptr<OverlappingNode<VertexT> > tree_iterator(m_root);
    tree_iterator->calculateOverlap(overlap_size);
}

template<typename VertexT>
void OverlappingKdTree<VertexT>::printKdTreeLeafs()
{
    boost::shared_ptr<OverlappingNode<VertexT> > tree_iterator(m_root);
    this->printKdTreeLeafs(tree_iterator,0);
}

template<typename VertexT>
void OverlappingKdTree<VertexT>::printKdTreeLeafs(boost::shared_ptr<OverlappingNode<VertexT> > node, unsigned int depth)
{
    if(node->isLeaf())
    {
        std::cout << "Node Size: " << node->getLeafSize() << std::endl;
        std::cout << "Overlap sizes: " ;
        for(unsigned int i=0; i<6; i++)
        {
            std::vector<VertexT> lap = node->getOverlap(i);
            std::cout << lap.size() << " " ;
        }
        std::cout << std::endl;
    } else {
        std::cout << "Depth: " << depth << " Split Left" << std::endl;
        this->printKdTreeLeafs(node->leftChild(), depth+1);

        std::cout << "Depth: " << depth << " Split Right" << std::endl;
        this->printKdTreeLeafs(node->rightChild(), depth+1);
    }
}

template<typename VertexT>
void OverlappingKdTree<VertexT>::collectLeafs(boost::shared_ptr<OverlappingNode<VertexT> > node)
{
    if(node->isLeaf())
    {
        m_leafs.push_back(node.get());
    }else{
        collectLeafs(node->leftChild());
        collectLeafs(node->rightChild());
    }
}

template<typename VertexT>
size_t OverlappingKdTree<VertexT>::getNumLeafs()
{
    return m_leafs.size();
}

}//namespace lvr