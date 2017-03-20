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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <rply.h>
#include <lvr/reconstruction/cuda/calcNormalsCuda.h>
#include <lvr/io/ModelFactory.hpp>
#include "Options.hpp"

using namespace lvr;

int main(int argc, char** argv){
	cuda_normals::Options opt(argc, argv);
    cout << opt << endl;
    
    ModelPtr model = ModelFactory::readModel(opt.inputFile());
    
	
	size_t num_points;
	
	floatArr points;
	
	if (model->m_pointCloud ) {
		
		points = model->m_pointCloud->getPointArray(num_points);
		
		std::cout << "Read " << num_points << " Points" << std::endl;
		
	} else if (model->m_mesh) {
		
		points = model->m_mesh->getVertexArray(num_points);
		
		std::cout << "Read " << num_points << " Vertices" << std::endl;
		
	}
	
	
	floatArr normals = floatArr(new float[ num_points * 3 ]);
	
	std::cout << "Constructing kd-tree..." << std::endl;
	CalcNormalsCuda calculator(points, num_points);
    std::cout << "Finished kd-tree construction." << std::endl;
	
	
	calculator.setK(opt.kd());
	
	
	if(opt.useRansac()){
		calculator.setMethod("RANSAC");
	}else{
		calculator.setMethod("PCA");
	}
	
	
	
	calculator.setFlippoint(opt.flipx(), opt.flipy(), opt.flipz());
	
	
    std::cout << "Start Normal Calculation..." << std::endl;
	calculator.start();
    
	
	calculator.getNormals(normals);
	std::cout << "Finished Normal Calculation. " << std::endl;
	
	
	
    ModelPtr out_model;
    
    if( model->m_pointCloud) {
		
		PointBufferPtr buffer = model->m_pointCloud;
		buffer->setPointNormalArray( normals, num_points);
		out_model = ModelPtr(new Model(buffer));
		
	} else if( model->m_mesh) {
		
		MeshBufferPtr buffer = model->m_mesh;
		buffer->setVertexNormalArray( normals, num_points);
		out_model = ModelPtr(new Model(buffer));
		
	}
    
    ModelFactory::saveModel(out_model, "normals.ply");
    
}
