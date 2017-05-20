/*
 * This file is part of cudaNormals.
 *
 * cudaNormals is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Foobar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cudaNormals.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * calcNormalsCuda.h
 *
 * @author Alexander Mock
 * @author Matthias Greshake
 */

#include <boost/filesystem.hpp>

#include <lvr/reconstruction/cuda/CUDANormals.hpp>
#include <lvr/io/ModelFactory.hpp>
#include <lvr/io/Timestamp.hpp>
#include <lvr/io/IOUtils.hpp>
#include "Options.hpp"


using namespace lvr;

void computeNormals(string filename, cuda_normals::Options& opt, PointBufferPtr& buffer)
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

    floatArr normals = floatArr(new float[ num_points * 3 ]);

    cout << timestamp << "Constructing kd-tree..." << endl;
    CalcNormalsCuda calculator(points, num_points);
    cout << timestamp << "Finished kd-tree construction." << endl;

    calculator.setK(opt.kd());
    if(opt.useRansac())
    {
        calculator.setMethod("RANSAC");
    } else
    {
        calculator.setMethod("PCA");
    }
    calculator.setFlippoint(opt.flipx(), opt.flipy(), opt.flipz());

    cout << timestamp << "Start Normal Calculation..." << endl;
    calculator.start();

    calculator.getNormals(normals);
    cout << timestamp << "Finished Normal Calculation. " << endl;

    buffer->setPointArray(points, num_points);
    buffer->setPointNormalArray(normals, num_points);
}

int main(int argc, char** argv){

    
    cuda_normals::Options opt(argc, argv);
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

                    computeNormals(p.string(), opt, buffer);
                    transformPointCloudAndAppend(buffer, p, all_points, all_normals);

                }
            }
        }

        writePointsAndNormals(all_points, all_normals, opt.outputFile());
    }
    else
    {
        PointBufferPtr buffer(new PointBuffer);
        computeNormals(opt.inputFile(), opt, buffer);
        ModelPtr out_model(new Model(buffer));
        ModelFactory::saveModel(out_model, opt.outputFile());
    }





}
