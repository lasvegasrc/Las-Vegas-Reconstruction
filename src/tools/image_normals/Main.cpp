/**
 * Copyright (C) 2017 Universität Osnabrück
 * This file is part of the LAS VEGAS Reconstruction Toolkit,
 *
 * LAS VEGAS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * LAS VEGAS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA
 */


// Program options for this tool

#include <lvr/io/IOUtils.hpp>
#include <lvr/reconstruction/ModelToImage.hpp>
#include <lvr/reconstruction/PanoramaNormals.hpp>
#include "Options.hpp"

#include <boost/filesystem.hpp>

#include <Eigen/Dense>

#include <fstream>
using std::ifstream;

using namespace lvr;


void computeNormals(string filename, image_normals::Options& opt, PointBufferPtr& buffer)
{
    ModelPtr model = ModelFactory::readModel(filename);

    // Determine coordinate system
    ModelToImage::CoordinateSystem system = ModelToImage::NATIVE;

    if(opt.coordinateSystem() == "SLAM6D")
    {
        system = ModelToImage::SLAM6D;
    }
    else if(opt.coordinateSystem() == "UOS")
    {
        system = ModelToImage::UOS;
    }

    ModelToImage mti(
                model->m_pointCloud,
                ModelToImage::CYLINDRICAL,
                opt.imageWidth(), opt.imageHeight(),
                opt.minZ(), opt.maxZ(),
                opt.minH(), opt.maxH(),
                opt.minV(), opt.maxV(),
                opt.optimize(), system);

    mti.writePGM(opt.imageFile(), 3000);

    PanoramaNormals normals(&mti);
    buffer = normals.computeNormals(opt.regionWidth(), opt.regionHeight(), false);
}



/**
 * @brief   Main entry point for the LSSR surface executable
 */
int main(int argc, char** argv)
{
    image_normals::Options opt(argc, argv);
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
                    PointBufferPtr buffer;
                    computeNormals(p.string(), opt, buffer);
                    transformPointCloudAndAppend(buffer, p, all_points, all_normals);

                }
            }
        }

        writePointsAndNormals(all_points, all_normals, opt.outputFile());
    }
    else
    {
        PointBufferPtr buffer;
        computeNormals(opt.inputFile(), opt, buffer);
        ModelPtr out_model(new Model(buffer));
        ModelFactory::saveModel(out_model, opt.outputFile());
    }

}

