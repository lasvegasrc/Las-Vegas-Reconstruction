#include <lvr/io/IOUtils.hpp>

namespace lvr
{

Eigen::Matrix4d buildTransformation(double* alignxf)
{
    Eigen::Matrix3d rotation;
    Eigen::Vector4d translation;

    rotation  << alignxf[0],  alignxf[4],  alignxf[8],
    alignxf[1],  alignxf[5],  alignxf[9],
    alignxf[2],  alignxf[6],  alignxf[10];

    translation << alignxf[12], alignxf[13], alignxf[14], 1.0;

    Eigen::Matrix4d transformation;
    transformation.setIdentity();
    transformation.block<3,3>(0,0) = rotation;
    transformation.rightCols<1>() = translation;

    return transformation;
}

Eigen::Matrix4d getTransformationFromPose(boost::filesystem::path& pose)
{
    std::ifstream poseIn(pose.c_str());
    if(poseIn.good())
    {
        double rPosTheta[3];
        double rPos[3];
        double alignxf[16];

        poseIn >> rPos[0] >> rPos[1] >> rPos[2];
        poseIn >> rPosTheta[0] >> rPosTheta[1] >> rPosTheta[2];

        rPosTheta[0] *= 0.0174533;
        rPosTheta[1] *= 0.0174533;
        rPosTheta[2] *= 0.0174533;

        double sx = sin(rPosTheta[0]);
        double cx = cos(rPosTheta[0]);
        double sy = sin(rPosTheta[1]);
        double cy = cos(rPosTheta[1]);
        double sz = sin(rPosTheta[2]);
        double cz = cos(rPosTheta[2]);

        alignxf[0]  = cy*cz;
        alignxf[1]  = sx*sy*cz + cx*sz;
        alignxf[2]  = -cx*sy*cz + sx*sz;
        alignxf[3]  = 0.0;
        alignxf[4]  = -cy*sz;
        alignxf[5]  = -sx*sy*sz + cx*cz;
        alignxf[6]  = cx*sy*sz + sx*cz;
        alignxf[7]  = 0.0;
        alignxf[8]  = sy;
        alignxf[9]  = -sx*cy;
        alignxf[10] = cx*cy;

        alignxf[11] = 0.0;

        alignxf[12] = rPos[0];
        alignxf[13] = rPos[1];
        alignxf[14] = rPos[2];
        alignxf[15] = 1;

        return buildTransformation(alignxf);
    }
    else
    {
        return Eigen::Matrix4d::Identity();
    }
}

Eigen::Matrix4d getTransformationFromFrames(boost::filesystem::path& frames)
{
    double alignxf[16];
    int color;

    std::ifstream in(frames.c_str());
    int c = 0;
    while(in.good())
    {
        c++;
        for(int i = 0; i < 16; i++)
        {
            in >> alignxf[i];
        }

        in >> color;

        if(!in.good())
        {
            c = 0;
            break;
        }
    }

    return buildTransformation(alignxf);
}


size_t countPointsInFile(boost::filesystem::path& inFile)
{
    std::ifstream in(inFile.c_str());
    std::cout << timestamp << "Counting points in " << inFile.filename().string() << "..." << std::endl;

    // Count lines in file
    size_t n_points = 0;
    char line[2048];
    while(in.good())
    {
        in.getline(line, 1024);
        n_points++;
    }
    in.close();

    std::cout << timestamp << "File " << inFile.filename().string() << " contains " << n_points << " points." << std::endl;

    return n_points;
}

void writeFrames(Eigen::Matrix4d transform, const boost::filesystem::path& framesOut)
{
    std::ofstream out(framesOut.c_str());

    // write the rotation matrix
    out << transform.col(0)(0) << " " << transform.col(0)(1) << " " << transform.col(0)(2) << " " << 0 << " "
        << transform.col(1)(0) << " " << transform.col(1)(1) << " " << transform.col(1)(2) << " " << 0 << " "
        << transform.col(2)(0) << " " << transform.col(2)(1) << " " << transform.col(2)(2) << " " << 0 << " ";

    // write the translation vector
    out << transform.col(3)(0) << " "
        << transform.col(3)(1) << " "
        << transform.col(3)(2) << " "
        << transform.col(3)(3);

    out.close();
}

size_t writeModel( ModelPtr model,const  boost::filesystem::path& outfile)
{
    size_t n_ip;
    floatArr arr = model->m_pointCloud->getPointArray(n_ip);

    ModelFactory::saveModel(model, outfile.string());

    return n_ip;
}

size_t writePointsToASCII(ModelPtr model, std::ofstream& out, bool nocolor)
{
    size_t n_ip, n_colors;

    floatArr arr = model->m_pointCloud->getPointArray(n_ip);

    ucharArr colors = model->m_pointCloud->getPointColorArray(n_colors);
    for(int a = 0; a < n_ip; a++)
    {
        out << arr[a * 3] << " " << arr[a * 3 + 1] << " " << arr[a * 3 + 2];

        if(n_colors && !(nocolor))
        {
            out << " " << (int)colors[a * 3] << " " << (int)colors[a * 3 + 1] << " " << (int)colors[a * 3 + 2];
        }
        out << std::endl;

    }

    return n_ip;
}

size_t getReductionFactor(ModelPtr model, size_t reduction)
{
    size_t n_points;
    floatArr arr = model->m_pointCloud->getPointArray(n_points);


    std::cout << timestamp << "Point cloud contains " << n_points << " points." << std::endl;

/*
     * If reduction is less than the number of points it will segfault
     * because the modulo operation is not defined for n mod 0
     * and we have to keep all points anyways.
     * Same if no targetSize was given.
     */
    if(reduction != 0)
    {
        if(reduction < n_points)
        {
            return (int)n_points / reduction;
        }
    }

    /* No reduction write all points */
    return 1;
}

size_t getReductionFactorASCII(boost::filesystem::path& inFile, size_t targetSize)
{
    /*
     * If reduction is less than the number of points it will segfault
     * because the modulo operation is not defined for n mod 0
     * and we have to keep all points anyways.
     * Same if no targetSize was given.
     */
    if(targetSize != 0)
    {
        // Count lines in file
        size_t n_points = countPointsInFile(inFile);

        if(targetSize < n_points)
        {
            return (int)n_points / targetSize;
        }
    }

    /* No reduction write all points */
    return 1;

}

void transformAndReducePointCloud(ModelPtr model, int modulo, int sx, int sy, int sz, int xPos, int yPos, int zPos)
{
    size_t n_ip, n_colors;
    size_t cntr = 0;

    floatArr arr = model->m_pointCloud->getPointArray(n_ip);
    ucharArr colors = model->m_pointCloud->getPointColorArray(n_colors);

    // Plus one because it might differ because of the 0-index
    // better waste memory for one float than having not enough space.
    // TO-DO think about exact calculation.
    size_t targetSize = (3 * ((n_ip)/modulo)) + modulo;
    floatArr points(new float[targetSize ]);
    ucharArr newColorsArr;

    if(n_colors)
    {
        newColorsArr = ucharArr(new unsigned char[targetSize]);
    }

    for(int i = 0; i < n_ip; i++)
    {
        if(i % modulo == 0)
        {
            if(sx != 1)
            {
                arr[i * 3] 		*= sx;
            }

            if(sy != 1)
            {
                arr[i * 3 + 1] 	*= sy;
            }

            if(sz != 1)
            {
                arr[i * 3 + 2] 	*= sz;
            }

            if((cntr * 3) < targetSize)
            {
                points[cntr * 3]     = arr[i * 3 + xPos];
                points[cntr * 3 + 1] = arr[i * 3 + yPos];
                points[cntr * 3 + 2] = arr[i * 3 + zPos];
            }
            else
            {
                std::cout << "The following is for debugging purpose" << std::endl;
                std::cout << "Cntr: " << (cntr * 3) << " targetSize: " << targetSize << std::endl;
                std::cout << "nip : " << n_ip << " modulo " << modulo << std::endl;
                break;
            }

            if(n_colors)
            {
                newColorsArr[cntr * 3]     = colors[i * 3];
                newColorsArr[cntr * 3 + 1] = colors[i * 3 + 1];
                newColorsArr[cntr * 3 + 2] = colors[i * 3 + 2];
            }

            cntr++;
        }
    }

    // Pass counter because it is the actual number of points used after reduction
    // it might be 1 less than the size
    model->m_pointCloud->setPointArray(points, cntr);

    if(n_colors)
    {
        model->m_pointCloud->setPointColorArray(newColorsArr, cntr);
    }
}


void transformPointCloud(ModelPtr model, Eigen::Matrix4d transformation)
{
    std::cout << timestamp << "Transforming model." << std::endl;
    size_t numPoints;

    floatArr arr = model->m_pointCloud->getPointArray(numPoints);

    for(int i = 0; i < numPoints; i++)
    {
        float x = arr[3 * i];
        float y = arr[3 * i + 1];
        float z = arr[3 * i + 2];

        Eigen::Vector4d v(x,y,z,1);
        Eigen::Vector4d tv = transformation * v;

        arr[3 * i]     = tv[0];
        arr[3 * i + 1] = tv[1];
        arr[3 * i + 2] = tv[2];
    }
}

} // namespace lvr
