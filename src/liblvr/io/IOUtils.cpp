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

} // namespace lvr
