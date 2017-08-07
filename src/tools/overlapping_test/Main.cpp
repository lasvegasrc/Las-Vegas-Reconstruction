#include <boost/filesystem.hpp>

#include <lvr/reconstruction/OverlappingKdTree.hpp>

#include <lvr/geometry/ColorVertex.hpp>
#include <lvr/io/ModelFactory.hpp>
#include <lvr/io/Timestamp.hpp>
#include <lvr/io/IOUtils.hpp>
#include "Options.hpp"


using namespace lvr;

typedef ColorVertex<double, unsigned char> cvertex;

void addPointFile(OverlappingKdTree<cvertex>& OlTree, string filename, overlapping_test::Options& opt, PointBufferPtr& buffer)
{
    ModelPtr model = ModelFactory::readModel(filename);
    size_t num_points;

    floatArr points;
    if (model && model->m_pointCloud )
    {
        points = model->m_pointCloud->getPointArray(num_points);
        cout << timestamp << "Read " << num_points << " points from " << filename << endl;
    }
    else
    {
        cout << timestamp << "Warning: No point cloud data found in " << filename << endl;
        return;
    }

    for(size_t i=0; i<num_points; i++)
    {
        cvertex in_point;
        in_point.x = points[i*3+0];
        in_point.y = points[i*3+1];
        in_point.z = points[i*3+2]; 
        OlTree.insert(in_point);
    }
    std::cout << "calculate Overlaps" << std::endl;
    double overlap_dist = 15;
    OlTree.calculateOverlaps(overlap_dist);
    OlTree.finishConstruction();

    OlTree.printKdTreeLeafs();

    std::cout << "Number of leafs: " << OlTree.getNumLeafs() << std::endl;
}

int main(int argc, char** argv){

    
    overlapping_test::Options opt(argc, argv);
    cout << opt << endl;

    
    unsigned int max_leaf_size = 200000;
    OverlappingKdTree<cvertex> OlTree(max_leaf_size);


    boost::filesystem::path inFile(opt.inputFile());

    if(boost::filesystem::is_directory(inFile))
    {
        vector<float> all_points;
        vector<float> all_normals;

        boost::filesystem::directory_iterator lastFile;
        for(boost::filesystem::directory_iterator it(inFile); it != lastFile; it++ )
        {
            boost::filesystem::path p = boost::filesystem::canonical(it->path());
            string currentFile = p.filename().string();

            if(string(p.extension().string().c_str()) == ".3d")
            {
                // Check for naming convention "scanxxx.3d"
                int num = 0;
                if(sscanf(currentFile.c_str(), "scan%3d", &num))
                {
                    cout << timestamp << "Processing " << p.string() << endl;
                    PointBufferPtr buffer(new PointBuffer);

                    addPointFile(OlTree, p.string(), opt, buffer);
                    transformPointCloudAndAppend(buffer, p, all_points, all_normals);

                }
            }
        }

        writePointsAndNormals(all_points, all_normals, opt.outputFile());
        
        

    }
    else
    {
        PointBufferPtr buffer(new PointBuffer);
        addPointFile(OlTree, opt.inputFile(), opt, buffer);

        ModelPtr out_model(new Model(buffer));
        ModelFactory::saveModel(out_model, opt.outputFile());
        
    }


}
