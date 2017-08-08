#include <boost/filesystem.hpp>

#include <lvr/reconstruction/OverlappingKdTree.hpp>

#include <lvr/geometry/ColorVertex.hpp>
#include <lvr/io/ModelFactory.hpp>
#include <lvr/io/Timestamp.hpp>
#include <lvr/io/IOUtils.hpp>
#include <fstream>
#include <string>
#include <algorithm>
#include "Options.hpp"


using namespace lvr;

typedef ColorVertex<double, unsigned char> cvertex;

void rgb(float minimum, float maximum, float value,
    unsigned int& r, unsigned int& g, unsigned int& b )
{
    float ratio = 2 * (value-minimum) / (maximum - minimum);
    b = int(std::max(float(0.0), 255*(1 - ratio)));
    r = int(std::max(float(0.0), 255*(ratio - 1)));
    g = 255 - b - r;
}
    

void convertToRgb(unsigned int minimum, unsigned int maximum, unsigned int value, 
        unsigned int& r, unsigned int& g, unsigned int& b )
{
    unsigned int halfmax = (minimum + maximum)/2;
    if( minimum <= value <= halfmax)
    {
        r = 0;
        g = 255/(halfmax - minimum) * (value - minimum);
        b = 255 - 255/(halfmax - minimum) * (value - minimum);
    }else if(halfmax < value <= maximum)
    {
        r = 255/(maximum - halfmax) * (value - halfmax);
        g = 255 + -255/(maximum - halfmax)  * (value - halfmax);
        b = 0;
    }
}

void writeBlob(std::vector<std::vector<cvertex>* >& points, unsigned int number=0)
{
    std::string number_str = std::to_string(number);

    if(number_str.length() == 1)
    {
        number_str = std::string("00") + number_str;
    }else if(number_str.length() == 2)
    {
        number_str = std::string("0") + number_str;
    }

    std::string filename = std::string("scan") + number_str + std::string(".blob");
    ofstream ofile(filename.c_str(), ios::out | ios::binary);

    unsigned int num_points = 0;
    for(unsigned int i=0; i< points.size(); i++)
    {
        num_points += points[i]->size();
    }

    ofile << num_points << std::endl;
    for(unsigned int i=0; i< points.size(); i++)
    {
        unsigned int num_points = points[i]->size();
        if(num_points > 0)
        {
            std::cout << "Writing " << num_points << " to " << filename << std::endl;
            ofile << num_points << std::endl;
            ofile.write( (char*)&((*points[i])[0]), num_points * sizeof(cvertex) );
        }
        
    }

    ofile.close();
}

void readBlob(unsigned int number, std::vector<cvertex>& out_points, unsigned int& num_inner_points )
{

    std::string number_str = std::to_string(number);

    if(number_str.length() == 1)
    {
        number_str = std::string("00") + number_str;
    }else if(number_str.length() == 2)
    {
        number_str = std::string("0") + number_str;
    }

    std::string filename = std::string("scan") + number_str + std::string(".blob");

    out_points.clear();
    ifstream infile(filename.c_str(), ios::in | ios::binary);
    unsigned int num_points_total;
    infile >> num_points_total;

    std::cout << "-- Read " << num_points_total << " Points in total" << std::endl;

    out_points.resize(num_points_total);



    char ch; // endl
    char point_size_buf[4];

    if(infile.is_open())
    {
        unsigned int i = 0;
        unsigned int index = 0;
        while(!infile.eof())
        {
            unsigned int point_size;
            infile >> point_size;

            if(i == 0)
            {
                num_inner_points = point_size;
            }

            infile.read((char*)&ch, 1);
            std::cout << "--- Read Blob " << i  << ", Num Points: " << point_size << std::endl;
            
            infile.read((char*)&(out_points[index]), sizeof(cvertex) * point_size );

            i++;
            index += point_size;
        }

        infile.close();
    }
    
}

void savePly(std::vector<std::vector<cvertex>* >& points, unsigned int number=0 )
{
    unsigned int num_points = 0;
    for(unsigned int i=0; i< points.size(); i++)
    {
        num_points += points[i]->size();
    }

    std::string filename = std::string("debug") + std::to_string(number) + std::string(".ply");

    std::cout << "Writing " << num_points << " to " << filename << std::endl;

    ofstream myfile;
    
    
    myfile.open(filename.c_str());
    myfile << "ply" << std::endl;
    myfile << "format ascii 1.0" << std::endl;
    myfile << "element vertex " << num_points << std::endl;
    myfile << "property float32 x" << std::endl;
    myfile << "property float32 y" << std::endl;
    myfile << "property float32 z" << std::endl;
    myfile << "property uchar red" << std::endl;                   
    myfile << "property uchar green" << std::endl;
    myfile << "property uchar blue" << std::endl;
    myfile << "end_header" << std::endl;

    for(unsigned int i=0; i<points.size(); i++)
    {
        for(unsigned int j=0; j<points[i]->size(); j++)
        {
            
            unsigned int r,g,b;
            //convertToRgb(0, points.size(), i, r, g, b);
            rgb(0, points.size(), i, r, g, b);

            myfile << (*points[i])[j][0] << " " 
                   << (*points[i])[j][1] << " " 
                   << (*points[i])[j][2] << " "
                   << r << " "
                   << g << " "
                   << b << std::endl;
        }
    }

    myfile.close();

}

unsigned int splitPointFile( string filename, overlapping_test::Options& opt, PointBufferPtr& buffer)
{
    unsigned int max_leaf_size = opt.maxLeafSize();
    OverlappingKdTree<cvertex> OlTree(max_leaf_size);

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
        return 0;
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
    double overlap_dist = opt.overlapSize();
    OlTree.calculateOverlaps(overlap_dist);
    OlTree.finishConstruction();

    OlTree.printKdTreeLeafs();

    std::cout << "Number of leafs: " << OlTree.getNumLeafs() << std::endl;

    std::vector<std::vector<std::vector< cvertex>* > > splitted_points;
    splitted_points = OlTree.getSplittedPoints();

    std::cout << std::endl;
    std::cout << "Splitted Points:" << std::endl;
    for(unsigned int i=0; i<splitted_points.size(); i++)
    {
        // debug
        writeBlob(splitted_points[i], i);
        
        std::cout << " " << i << std::endl;
        for(unsigned int j=0; j<splitted_points[i].size(); j++)
        {
            std::cout << "  " << j << ": " << splitted_points[i][j]->size() << std::endl;
        }
    }

    return splitted_points.size();
}

int main(int argc, char** argv){

    
    overlapping_test::Options opt(argc, argv);
    cout << opt << endl;

    
    


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

                    splitPointFile( p.string(), opt, buffer);
                    transformPointCloudAndAppend(buffer, p, all_points, all_normals);

                }
            }
        }

        writePointsAndNormals(all_points, all_normals, opt.outputFile());
        
    }
    else
    {
        PointBufferPtr buffer(new PointBuffer);
        unsigned int num_files = splitPointFile( opt.inputFile(), opt, buffer);


        for(unsigned int i = 0; i< num_files; i++)
        {
            std::cout << "Read File " << i << std::endl;
            std::vector< cvertex> points;
            unsigned int num_inner_points;
            readBlob(i, points, num_inner_points);
            std::cout << "Number of Points: " << points.size() << std::endl;
            std::cout << "   Number of inner Points: " << num_inner_points << std::endl;
        }
        

        ModelPtr out_model(new Model(buffer));
        ModelFactory::saveModel(out_model, opt.outputFile());
        
    }


}
