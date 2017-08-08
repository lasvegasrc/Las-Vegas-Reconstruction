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

    this->collectNeighbours();
}

template<typename VertexT>
void OverlappingKdTree<VertexT>::getNodeOverlapsWithId(
    boost::shared_ptr<OverlappingNode<VertexT> > node,
    unsigned int lap_id,
    std::vector< std::vector<VertexT>* >& output)
{
    if(node->isLeaf())
    {
        output.push_back( node->getOverlap(lap_id) );
    }else{
        getNodeOverlapsWithId(node->m_left, lap_id, output);
        getNodeOverlapsWithId(node->m_right, lap_id, output);
    }
}

template<typename VertexT>
void OverlappingKdTree<VertexT>::collectNeighbours()
{
    for(auto it = m_leafs.begin(); it != m_leafs.end(); ++it)
    {
        std::vector<std::vector<VertexT>* > points;
        points.push_back((*it)->m_points.get());
        // and neighbours overlaps
        
        for(auto map_it = (*it)->m_neighbours.begin();
                map_it != (*it)->m_neighbours.end();
                map_it++ )
        {
            std::vector< std::vector<VertexT>* > neighbours_overlaps;
            unsigned int dim = map_it->first/2;
            unsigned int side = (map_it->first+1)%2;
            getNodeOverlapsWithId(map_it->second, dim*2+side, neighbours_overlaps);
            for(int i = 0; i<neighbours_overlaps.size(); i++)
            {
                points.push_back(neighbours_overlaps[i]);
            }
        }
        m_points.push_back(points);
    }
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
void OverlappingKdTree<VertexT>::printNeighbours(boost::shared_ptr<OverlappingNode<VertexT> > node)
{
    if(node->isLeaf())
    {
        std::cout << node->m_points->size() << "|"; 
    }else{
        printNeighbours(node->m_left);
        printNeighbours(node->m_right);
    }
}

template<typename VertexT>
void OverlappingKdTree<VertexT>::printKdTreeLeafs(boost::shared_ptr<OverlappingNode<VertexT> > node, unsigned int depth)
{
    if(node->isLeaf())
    {
        std::cout << "---Node Size: " << node->getLeafSize() << std::endl;
        std::cout << "---Overlap sizes: " ;
        for(unsigned int i=0; i<6; i++)
        {
            std::vector<VertexT>* lap = node->getOverlap(i);
            std::cout << lap->size() << " " ;
        }
        std::cout << std::endl;
        std::cout << "---Neighbor sizes: " ;
        for(auto it = node->m_neighbours.begin(); it != node->m_neighbours.end(); ++it)
        {
            if(it->second->isLeaf() )
            {
                std::cout << it->first << ": " << it->second->m_points->size() << ", ";
            }else{

                std::cout << it->first << ": ";
                printNeighbours(it->second);
            }
        }
        std::cout << std::endl;

        std::cout << "---Bounding Box: " << node->m_bounding_box << std::endl;
    } else {
        std::cout << "Depth: " << depth << ", Dim: "<< node->getSplitDimension()<< " Split Left" << std::endl;
        this->printKdTreeLeafs(node->leftChild(), depth+1);

        std::cout << "Depth: " << depth << ", Dim: "<< node->getSplitDimension() <<" Split Right" << std::endl;
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

template<typename VertexT>
std::vector< std::vector<std::vector< VertexT>* > >& OverlappingKdTree<VertexT>::getSplittedPoints()
{
    return m_points;
}


}//namespace lvr