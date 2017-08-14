#include <boost/filesystem.hpp>
#include <boost/shared_array.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_lit.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include <lvr/reconstruction/OverlappingKdTree.hpp>

#include <lvr/geometry/ColorVertex.hpp>
#include <lvr/io/ModelFactory.hpp>
#include <lvr/io/Progress.hpp>
#include <lvr/io/Timestamp.hpp>
#include <lvr/io/IOUtils.hpp>

#include <lvr/reconstruction/cuda/CudaSurface.hpp>

#include <fstream>
#include <string>
#include <algorithm>
#include "Options.hpp"


using namespace lvr;

typedef ColorVertex<double, unsigned char> cvertex;

void rgb(float minimum, float maximum, float value,
    size_t& r, size_t& g, size_t& b )
{
    float ratio = 2 * (value-minimum) / (maximum - minimum);
    b = int(std::max(float(0.0), 255*(1 - ratio)));
    r = int(std::max(float(0.0), 255*(ratio - 1)));
    g = 255 - b - r;
}
    

void convertToRgb(size_t minimum, size_t maximum, size_t value,
        size_t& r, size_t& g, size_t& b )
{
    size_t halfmax = (minimum + maximum)/2;
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

void writeBlob(std::vector<std::vector<cvertex>* >& points, size_t number=0)
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

    size_t num_points = 0;
    for(size_t i=0; i< points.size(); i++)
    {
        num_points += points[i]->size();
    }

    ofile << num_points << std::endl;
    for(size_t i=0; i< points.size(); i++)
    {
        size_t num_points = points[i]->size();
        if(num_points > 0)
        {
            std::cout << "Writing " << num_points << " to " << filename << std::endl;
            ofile << num_points << std::endl;
            ofile.write( (char*)&((*points[i])[0]), num_points * sizeof(cvertex) );
        }
        
    }

    ofile.close();
}

void readBlob(size_t number, std::vector<cvertex>& out_points, size_t& num_inner_points )
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
    size_t num_points_total;
    infile >> num_points_total;

    std::cout << "-- Read " << num_points_total << " Points in total" << std::endl;

    out_points.resize(num_points_total);



    char ch; // endl
    char point_size_buf[4];

    if(infile.is_open())
    {
        size_t i = 0;
        size_t index = 0;
        while(!infile.eof())
        {
            size_t point_size;
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

void savePly(std::vector<std::vector<cvertex>* >& points, size_t number=0 )
{
    size_t num_points = 0;
    for(size_t i=0; i< points.size(); i++)
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

    for(size_t i=0; i<points.size(); i++)
    {
        for(size_t j=0; j<points[i]->size(); j++)
        {
            
            size_t r,g,b;
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

void addToTree(OverlappingKdTree<cvertex>& tree,  boost::filesystem::path& transformFile)
{
    cout << timestamp << "Adding " << transformFile.c_str() << " to tree." << endl;

    char frames[2048];
    char pose[2014];

    sprintf(frames, "%s/%s.frames", transformFile.parent_path().c_str(), transformFile.stem().c_str());
    sprintf(pose, "%s/%s.pose", transformFile.parent_path().c_str(), transformFile.stem().c_str());

    boost::filesystem::path framesPath(frames);
    boost::filesystem::path posePath(pose);


    Eigen::Matrix4d transform = Eigen::Matrix4d::Identity();

    if(boost::filesystem::exists(framesPath))
    {
       std::cout << timestamp << "Transforming according to " << framesPath.filename() << std::endl;
       transform = getTransformationFromFrames(framesPath);
    }
    else if(boost::filesystem::exists(posePath))
    {
       std::cout << timestamp << "Transforming according to " << posePath.filename() << std::endl;
       transform = getTransformationFromFrames(posePath);
    }
    else
    {
       std::cout << timestamp << "Warning: found no transformation for " << transformFile.filename() << std::endl;
    }

    cout << timestamp << "Reading point cloud" << endl;

    ModelPtr model = ModelFactory::readModel(transformFile.c_str());
    PointBufferPtr buffer = model->m_pointCloud;

    size_t n_points;
    size_t n_colors;

    floatArr points = buffer->getPointArray(n_points);
    ucharArr colors = buffer->getPointColorArray(n_colors);

    string comment = timestamp.getElapsedTime() + "Transforming points ";
    ProgressBar progress(n_points, comment);


    for(size_t i = 0; i < n_points; i++)
    {
        float x = points[3 * i];
        float y = points[3 * i + 1];
        float z = points[3 * i + 2];

        Eigen::Vector4d v(x,y,z,1);
        Eigen::Vector4d tv = transform * v;

        cvertex cv;
        cv.x = tv[0];
        cv.y = tv[1];
        cv.z = tv[2];

        if(n_colors)
        {
            cv.r = colors[3 * i];
            cv.g = colors[3 * i + 1];
            cv.b = colors[3 * i + 2];
        }
        else
        {
            cv.r = 0;
            cv.g = 255;
            cv.b = 0;
        }

        tree.insert(cv);
        ++progress;
    }
    cout << endl;

}

void getPointBufferFromBlob(PointBufferPtr& buffer, int blob_number, size_t& num_inner_points)
{
    char n[1025];
    sprintf(n, "scan%03d.blob", blob_number);

    num_inner_points = 0;

    cout << timestamp << "Parsing " << n << endl;
//    ifstream infile(n, ios::in | ios::binary);

//    size_t num_points_total;
//    infile >> num_points_total;


//    char ch; // endl
//    char point_size_buf[4];

//    if(infile.is_open())
//    {
//        size_t i = 0;
//        size_t index = 0;

//        while(!infile.eof())
//        {
//            size_t point_size;
//            infile >> point_size;

//            // Check needed. Last read could have invalidated the stream
//            if(!infile.eof())
//            {
//                cvertex tmp[10];

//                if(i == 0)
//                {
//                    num_inner_points = point_size;
//                }

//                infile.read((char*)&ch, 1);
//                std::cout << "--- Read Blob " << i  << ", Num Points: " << point_size << std::endl;

//                int n_chunks = point_size / 10;
//                int n_rest = point_size % 10;

//                for(int a = 0; a < n_chunks; a++)
//                {
//                    infile.read((char*)tmp, sizeof(cvertex) * 10 );
//                }

//                for(int b = 0; b < n_rest; b++)
//                {
//                    infile.read((char*)tmp, sizeof(cvertex) * n_rest );
//                }


//                //infile.read((char*)&(out_points[index]), sizeof(cvertex) * point_size );

//                i++;

//                cout << sizeof(cvertex) * point_size << " " << infile.gcount() << endl;
//                index += point_size;
//            }

//        }
//        infile.sync();
//        infile.close();
//    }
    size_t num_points_total = 0;
    char newline;
    FILE* f = fopen(n, "r");
    if(f)
    {
        fscanf(f, "%zu", &num_points_total);
        cout << "Num points total: " << num_points_total << endl;
        cvertex* out_points = new cvertex[num_points_total];
        cvertex  pt;

        cout << sizeof(cvertex) << endl;

        if(num_points_total)
        {
            size_t i = 0;
            size_t index = 0;
            size_t pos = 0;
            do
            {
                size_t point_size = 0;
                fscanf(f, "%zu", &point_size);

                if(feof(f)) break;
                fread(&newline, 1, 1, f);

                if(i == 0)
                {
                    num_inner_points = point_size;
                }

                std::cout << "--- Read Blob " << i  << ", Num Points: " << point_size << std::endl;

                for(int g = 0; g < point_size; g++)
                {
                    fread(&pt, sizeof(cvertex), 1, f);
                    out_points[pos] = pt;
                    pos++;
                }

                i++;
                index+= point_size;
            }
            while (!feof(f));

            cout << timestamp << "Converting to point buffer" << endl;

            floatArr b_points(new float[3 * num_points_total]);
            ucharArr b_colors(new unsigned char[3 * num_points_total]);

            for(size_t j = 0; j < num_points_total; j++)
            {
                b_points[3 * j]     = out_points[j].x;
                b_points[3 * j + 1] = out_points[j].y;
                b_points[3 * j + 2] = out_points[j].z;

                b_colors[3 * j]     = out_points[j].r;
                b_colors[3 * j + 1] = out_points[j].g;
                b_colors[3 * j + 2] = out_points[j].b;

            }

            buffer->setPointArray(b_points, num_points_total);
            buffer->setPointColorArray(b_colors, num_points_total);
        }

        //cout << out_points << " " << ptr << endl;
        delete[] out_points;
    }
}

void calculateNormalsCUDA(size_t n_blobs, size_t kn, size_t ki)
{
    size_t n_all = 0;
    for(int i = 0; i < n_blobs; i++)
    {
        size_t n_inner = 0;

        PointBufferPtr buffer(new PointBuffer);
        getPointBufferFromBlob(buffer, i, n_inner);

        size_t n_points;
        floatArr points = buffer->getPointArray(n_points);

        size_t n_colors;
        ucharArr colors = buffer->getPointColorArray(n_colors);

        floatArr normals = floatArr(new float[n_points * 3]);


        cout << timestamp << "Generating CUDA surface" << endl;
        CudaSurface gpu_surface(points, n_points);
        gpu_surface.setKn(kn);
        gpu_surface.setKi(ki);
        gpu_surface.setMethod("PCA");

        cout << timestamp << "Starting normal calculation" << endl;
        gpu_surface.calculateNormals();

        gpu_surface.getNormals(normals);
        cout << timestamp << "Finished Normal Calculation. " << endl;

        gpu_surface.freeGPU();

        cout << timestamp << "Generating tmp data for " << n_inner << " inner points" << endl;

        size_t size = (6 * sizeof(float) + 3) * n_inner;

        unsigned char* arr = new unsigned char[size];
        unsigned char* p = &arr[0];

        n_all += n_inner;

        for(size_t j = 0; j < n_inner; j++)
        {
            float* t = (float*)p;
            *t = points[3 * j];
            t++;
            *t = points[3 * j + 1];
            t++;
            *t = points[3 * j + 2];
            t++;

            unsigned char* u = (unsigned char*)t;
            *u = colors[3 * j];
            u++;
            *u = colors[3 * j + 1];
            u++;
            *u = colors[3 * j + 2];
            u++;

            t = (float*)u;
            *t = normals[3 * j];
            t++;
            *t = normals[3 * j + 1];
            t++;
            *t = normals[3 * j + 2];
            t++;

            p = (unsigned char*)t;
        }

        char f_buf[1024];
        sprintf(f_buf,"chunk%03d.tmp", i);
        cout << "Writing tmp data to " << f_buf << endl;

        ofstream out(f_buf, ios::out | ios::binary);
        out.write( (char*)(arr), size);
        out.close();

        delete[] arr;
        buffer.reset();
        normals.reset();
        points.reset();
        colors.reset();
    }

    cout << timestamp << "Generating PLY Header" << endl;

    ofstream ply_out("normals.ply", ios::out | ios::binary);
    ply_out << "ply" << endl;
    ply_out << "format binary_little_endian 1.0" << std::endl;
    ply_out << "element vertex " << n_all << std::endl;
    ply_out << "property float x" << std::endl;
    ply_out << "property float y" << std::endl;
    ply_out << "property float z" << std::endl;
    ply_out << "property uchar red" << std::endl;
    ply_out << "property uchar green" << std::endl;
    ply_out << "property uchar blue" << std::endl;
    ply_out << "property float nx" << std::endl;
    ply_out << "property float ny" << std::endl;
    ply_out << "property float nz" << std::endl;
    ply_out << "end_header" << std::endl;

    cout << timestamp << "Concatenating chunks" << endl;

    for(int i = 0; i < n_blobs; i++)
    {
        char f_buf[1024];
        sprintf(f_buf,"chunk%03d.tmp", i);
        ifstream in(f_buf, ios::in | ios::binary);
        cout << timestamp << "Adding " << f_buf << endl;
        ply_out << in.rdbuf();
        in.close();
    }
    ply_out.close();
}

size_t splitPointFile( string filename, overlapping_test::Options& opt, PointBufferPtr& buffer)
{
    size_t max_leaf_size = opt.maxLeafSize();
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
    for(size_t i=0; i<splitted_points.size(); i++)
    {
        // debug
        writeBlob(splitted_points[i], i);
        
        std::cout << " " << i << std::endl;
        for(size_t j=0; j<splitted_points[i].size(); j++)
        {
            std::cout << "  " << j << ": " << splitted_points[i][j]->size() << std::endl;
        }
    }

    return splitted_points.size();
}

template <typename Iterator>
bool parse_filename(Iterator first, Iterator last, int& i)
{

using boost::spirit::qi::lit;
using boost::spirit::qi::uint_parser;
using boost::spirit::qi::parse;
using boost::spirit::qi::_1;
using boost::phoenix::ref;

uint_parser<unsigned, 10, 3, 3> uint_3_d;

bool r = parse(
        first,                          /*< start iterator >*/
        last,                           /*< end iterator >*/
        ((lit("scan")|lit("Scan")) >> uint_3_d[ref(i) = _1])   /*< the parser >*/
        );

if (first != last) // fail if we did not get a full match
    return false;
return r;
}

bool sortPaths(boost::filesystem::path firstScan, boost::filesystem::path secScan)
{
    std::string firstStem = firstScan.stem().string();
    std::string secStem   = secScan.stem().string();

    int i = 0;
    int j = 0;

    bool first = parse_filename(firstStem.begin(), firstStem.end(), i);
    bool sec = parse_filename(secStem.begin(), secStem.end(), j);

    if(first && sec)
    {
        return (i < j);
    }
    else
    {
        // this causes non valid files being at the beginning of the vector.
        if(sec)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


int main(int argc, char** argv){

    
    overlapping_test::Options opt(argc, argv);
    cout << opt << endl;

    if(opt.blobMode())
    {
        size_t n = opt.nBlobs();
        calculateNormalsCUDA(n, opt.kn(), opt.ki());
    }
    else
    {

        boost::filesystem::path inFile(opt.inputFile());

        if(boost::filesystem::is_directory(inFile))
        {
            size_t max_leaf_size = opt.maxLeafSize();
            OverlappingKdTree<cvertex> olTree(max_leaf_size);

            vector<boost::filesystem::path> files;
            boost::filesystem::directory_iterator lastFile;
            for(boost::filesystem::directory_iterator it(inFile); it != lastFile; it++ )
            {
                boost::filesystem::path p = boost::filesystem::canonical(it->path());
                string currentFile = p.filename().string();

                if(string(p.extension().string().c_str()) == ".3d")
                {
                    files.push_back(p);
                    //addToTree(olTree, p);
                }


            }

            // Sort files and determine first and last file to parse
            std::sort(files.begin(), files.end(), sortPaths);

            size_t start = 0;
            size_t end = files.size();

            if(opt.first() > start && opt.first() < end)
            {
                start = opt.first();
            }

            if(opt.last() < end && opt.last() > start)
            {
                end = opt.last();
            }

            cout << start << " " << end << endl;

            for(size_t i = start; i < end; i++)
            {
                boost::filesystem::path p = files[i];
                addToTree(olTree, p);
            }


            double overlap_dist = opt.overlapSize();
            cout << "Finishing tree with overlap distance " << overlap_dist << endl;

            olTree.calculateOverlaps(overlap_dist);
            olTree.finishConstruction();

            olTree.printKdTreeLeafs();

            cout << timestamp << " Writing blobs ..." << endl;
            std::vector<std::vector<std::vector< cvertex>* > > splitted_points;
            splitted_points = olTree.getSplittedPoints();

            std::cout << std::endl;
            std::cout << "Splitted Points:" << std::endl;
            for(size_t i=0; i<splitted_points.size(); i++)
            {
                // debug
                writeBlob(splitted_points[i], i);

                std::cout << " " << i << std::endl;
                for(size_t j=0; j<splitted_points[i].size(); j++)
                {
                    std::cout << "  " << j << ": " << splitted_points[i][j]->size() << std::endl;
                }
            }
            cout << timestamp << "done..." << endl;
        }
        else
        {
            PointBufferPtr buffer(new PointBuffer);
            size_t num_files = splitPointFile( opt.inputFile(), opt, buffer);


            for(size_t i = 0; i< num_files; i++)
            {
                std::cout << "Read File " << i << std::endl;
                std::vector< cvertex> points;
                size_t num_inner_points;
                readBlob(i, points, num_inner_points);
                std::cout << "Number of Points: " << points.size() << std::endl;
                std::cout << "   Number of inner Points: " << num_inner_points << std::endl;
            }


            ModelPtr out_model(new Model(buffer));
            ModelFactory::saveModel(out_model, opt.outputFile());

        }

    }
}
