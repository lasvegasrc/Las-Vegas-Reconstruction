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

typedef boost::shared_array<unsigned int> uintArr;


typedef boost::shared_array<float> floatArr;


typedef boost::shared_array<unsigned char> ucharArr;


int readVertexCb( p_ply_argument argument )
{
    float ** ptr;
    ply_get_argument_user_data( argument, (void **) &ptr, NULL );
    **ptr = ply_get_argument_value( argument );
    (*ptr)++;
    return 1;

}

int readColorCb( p_ply_argument argument )
{

    uint8_t ** color;
    ply_get_argument_user_data( argument, (void **) &color, NULL );
    **color = ply_get_argument_value( argument );
    (*color)++;
    return 1;

}

int readFaceCb( p_ply_argument argument )
{

    unsigned int ** face;
    long int length, value_index;
    ply_get_argument_user_data( argument, (void **) &face, NULL );
    ply_get_argument_property( argument, NULL, &length, &value_index );
    if ( value_index < 0 )
    {
        /* We got info about amount of face vertices. */
        if ( ply_get_argument_value( argument ) == 3 )
        {
            return 1;
        }
        std::cerr  << "Mesh is not a triangle mesh." << std::endl;
        return 0;
    }
    **face = ply_get_argument_value( argument );
    (*face)++;

    return 1;

}

void readPlyFile(PointArray& V, size_t& num_points, const char* filename ,float scale=1.0, int number_of_data = -1,bool readColor =true, bool readConfidence=true,
    bool readIntensity=true, bool readNormals=true, bool readFaces=true){
	/* Start reading new PLY */
    p_ply ply = ply_open( filename, NULL, 0, NULL );

    if ( !ply )
    {
       std::cerr  << "Could not open »" << filename << "«."
           << std::endl;
        return ;
    }
    if ( !ply_read_header( ply ) )
    {
       std::cerr  << "Could not read header." << std::endl;
        return;
    }

    /* Check if there are vertices and get the amount of vertices. */
    char buf[256] = "";
    const char * name = buf;
    long int n;
    p_ply_element elem  = NULL;

    // Buffer count variables
    size_t numVertices              = 0;
    size_t numVertexColors          = 0;
    size_t numVertexConfidences     = 0;
    size_t numVertexIntensities     = 0;
    size_t numVertexNormals         = 0;

    size_t numPoints                = 0;
    size_t numPointColors           = 0;
    size_t numPointConfidence       = 0;
    size_t numPointIntensities      = 0;
    size_t numPointNormals          = 0;
    size_t numFaces                 = 0;

    while ( ( elem = ply_get_next_element( ply, elem ) ) )
    {
        ply_get_element_info( elem, &name, &n );
        if ( !strcmp( name, "vertex" ) )
        {
            numVertices = n;
            p_ply_property prop = NULL;
            while ( ( prop = ply_get_next_property( elem, prop ) ) )
            {
                ply_get_property_info( prop, &name, NULL, NULL, NULL );
                if ( !strcmp( name, "red" ) && readColor )
                {
                    /* We have color information */
                    numVertexColors = n;
                }
                else if ( !strcmp( name, "confidence" ) && readConfidence )
                {
                    /* We have confidence information */
                    numVertexConfidences = n;
                }
                else if ( !strcmp( name, "intensity" ) && readIntensity )
                {
                    /* We have intensity information */
                    numVertexIntensities = n;
                }
                else if ( !strcmp( name, "nx" ) && readNormals )
                {
                    /* We have normals */
                    numVertexNormals = n;
                }
            }
        }
        else if ( !strcmp( name, "point" ) )
        {
            numPoints = n;
            p_ply_property prop = NULL;
            while ( ( prop = ply_get_next_property( elem, prop ) ) )
            {
                ply_get_property_info( prop, &name, NULL, NULL, NULL );
                if ( !strcmp( name, "red" ) && readColor )
                {
                    /* We have color information */
                    numPointColors = n;
                }
                else if ( !strcmp( name, "confidence" ) && readConfidence )
                {
                    /* We have confidence information */
                    numPointConfidence = n;
                }
                else if ( !strcmp( name, "intensity" ) && readIntensity )
                {
                    /* We have intensity information */
                    numPointIntensities = n;
                }
                else if ( !strcmp( name, "nx" ) && readNormals )
                {
                    /* We have normals */
                    numPointNormals = n;
                }
            }
        }
        else if ( !strcmp( name, "face" ) && readFaces )
        {
            numFaces = n;
        }
    }

    if ( !( numVertices || numPoints ) )
    {
        std::cout << "Neither vertices nor points in ply."
            << std::endl;
        return ;
    }

    // Buffers
    floatArr vertices;
    floatArr vertexConfidence;
    floatArr vertexIntensity;
    floatArr vertexNormals;
    floatArr points;
    floatArr pointConfidences;
    floatArr pointIntensities;
    floatArr pointNormals;

    ucharArr pointColors;
	ucharArr vertexColors;
    uintArr  faceIndices;


    /* Allocate memory. */
    if ( numVertices )
    {
        vertices = floatArr( new float[ numVertices * 3 ] );
    }
    if ( numVertexColors )
    {
        vertexColors = ucharArr( new unsigned char[ numVertices * 3 ] );
    }
    if ( numVertexConfidences )
    {
        vertexConfidence = floatArr( new float[ numVertices ] );
    }
    if ( numVertexIntensities )
    {
        vertexIntensity = floatArr( new float[ numVertices ] );
    }
    if ( numVertexNormals )
    {
        vertexNormals = floatArr( new float[ 3 * numVertices ] );
    }
    if ( numFaces )
    {
        faceIndices = uintArr( new unsigned int[ numFaces * 3 ] );
    }
    if ( numPoints )
    {
        points = floatArr( new float[ numPoints * 3 ] );
    }
    if ( numPointColors )
    {
        pointColors = ucharArr( new unsigned char[ numPoints * 3 ] );
    }
    if ( numPointConfidence )
    {
        pointConfidences = floatArr( new float[numPoints] );
    }
    if ( numPointIntensities )
    {
        pointIntensities = floatArr( new float[numPoints] );
    }
    if ( numPointNormals )
    {
        pointNormals = floatArr( new float[ 3 * numPoints ] );
    }


    float*        vertex            = vertices.get();
    uint8_t* 	  vertex_color      = vertexColors.get();
    float*        vertex_confidence = vertexConfidence.get();
    float*        vertex_intensity  = vertexIntensity.get();
    float*        vertex_normal     = vertexNormals.get();
    unsigned int* face              = faceIndices.get();
    float*        point             = points.get();
    uint8_t*      point_color       = pointColors.get();
    float*        point_confidence  = pointConfidences.get();
    float*        point_intensity   = pointIntensities.get();
    float*        point_normal      = pointNormals.get();


    /* Set callbacks. */
    if ( vertex )
    {
        ply_set_read_cb( ply, "vertex", "x", readVertexCb, &vertex, 0 );
        ply_set_read_cb( ply, "vertex", "y", readVertexCb, &vertex, 0 );
        ply_set_read_cb( ply, "vertex", "z", readVertexCb, &vertex, 1 );
    }
    if ( vertex_color )
    {
        ply_set_read_cb( ply, "vertex", "red",   readColorCb,  &vertex_color,  0 );
        ply_set_read_cb( ply, "vertex", "green", readColorCb,  &vertex_color,  0 );
        ply_set_read_cb( ply, "vertex", "blue",  readColorCb,  &vertex_color,  1 );
    }
    if ( vertex_confidence )
    {
        ply_set_read_cb( ply, "vertex", "confidence", readVertexCb, &vertex_confidence, 1 );
    }
    if ( vertex_intensity )
    {
        ply_set_read_cb( ply, "vertex", "intensity", readVertexCb, &vertex_intensity, 1 );
    }
    if ( vertex_normal )
    {
        ply_set_read_cb( ply, "vertex", "nx", readVertexCb, &vertex_normal, 0 );
        ply_set_read_cb( ply, "vertex", "ny", readVertexCb, &vertex_normal, 0 );
        ply_set_read_cb( ply, "vertex", "nz", readVertexCb, &vertex_normal, 1 );
    }

    if ( face )
    {
        ply_set_read_cb( ply, "face", "vertex_indices", readFaceCb, &face, 0 );
        ply_set_read_cb( ply, "face", "vertex_index", readFaceCb, &face, 0 );
    }

    if ( point )
    {
        ply_set_read_cb( ply, "point", "x", readVertexCb, &point, 0 );
        ply_set_read_cb( ply, "point", "y", readVertexCb, &point, 0 );
        ply_set_read_cb( ply, "point", "z", readVertexCb, &point, 1 );
    }
    if ( point_color )
    {
        ply_set_read_cb( ply, "point", "red",   readColorCb,  &point_color,  0 );
        ply_set_read_cb( ply, "point", "green", readColorCb,  &point_color,  0 );
        ply_set_read_cb( ply, "point", "blue",  readColorCb,  &point_color,  1 );
    }
    if ( point_confidence )
    {
        ply_set_read_cb( ply, "point", "confidence", readVertexCb, &point_confidence, 1 );
    }
    if ( point_intensity )
    {
        ply_set_read_cb( ply, "point", "intensity", readVertexCb, &point_intensity, 1 );
    }
    if ( point_normal )
    {
        ply_set_read_cb( ply, "point", "nx", readVertexCb, &point_normal, 0 );
        ply_set_read_cb( ply, "point", "ny", readVertexCb, &point_normal, 0 );
        ply_set_read_cb( ply, "point", "nz", readVertexCb, &point_normal, 1 );
    }

    /* Read ply file. */
    if ( !ply_read( ply ) )
    {
        std::cerr << "Could not read »" << filename << "«."
            << std::endl;
    }



    ply_close( ply );
    
    
    if(vertices)
    {
		if(number_of_data > 0)
		{
			num_points = number_of_data;
			
			V.width = num_points;
			V.elements = (float*)malloc(number_of_data * 3 * sizeof(float) );
			
			for(int i=0; i<number_of_data; i++){
				int index = (int)((float)rand() * ((float)numVertices-1.0)/RAND_MAX   ) ;
				
				float x = vertices[index*3+0] * scale;
				float y = vertices[index*3+1] * scale;
				float z = vertices[index*3+2] * scale;
				
				//~ printf("%f %f %f \n", x,y,z);
				
				V.elements[i * 3] = x;
				V.elements[i * 3 + 1] = y;
				V.elements[i * 3 + 2] = z;
				
				
			}
		}else{
			
			V.width = numVertices;
			V.elements = (float*)malloc(numVertices*3*sizeof(float) );
			num_points = numVertices;
			for(int i=0; i<numVertices; i++){
			
				V.elements[i*3] = vertices[i*3+0]*scale;
				V.elements[i*3+1] = vertices[i*3+1]*scale;
				V.elements[i*3+2] = vertices[i*3+2]*scale;
			}
		}
    }else if(points) {
		
		if(number_of_data > 0)
		{
			num_points = number_of_data;
			V.elements = (float*)malloc(number_of_data * 3 * sizeof(float));
			
			for(int i=0; i<number_of_data; i++)
			{
				int index = (int)((float)rand() * ((float)numPoints-1.0) / RAND_MAX   ) ;
				
				float x = points[index*3+0] * scale;
				float y = points[index*3+1] * scale;
				float z = points[index*3+2] * scale;
				
				//~ printf("%f %f %f \n", x,y,z);
				
				V.elements[i*3] = x;
				V.elements[i*3+1] = y;
				V.elements[i*3+2] = z;
				
			}
		}else{
			num_points = numPoints;
			V.elements = (float*)malloc(numPoints*3*sizeof(float) );
			
			for(int i=0; i<numPoints; i++)
			{
			
				V.elements[i*3] = points[i*3+0]*scale;
				V.elements[i*3+1] = points[i*3+1]*scale;
				V.elements[i*3+2] = points[i*3+2]*scale;
			}
		}
	}
    

    
}

void writePlyFile(float* V, size_t m_numVertices, float* Result_Normals, size_t m_numVertexNormals, const char* filename, float scale=1.0){
	
	/* Handle options. */
    e_ply_storage_mode mode( PLY_LITTLE_ENDIAN );
    
    p_ply oply = ply_create( filename, mode, NULL, 0, NULL );
    if ( !oply )
    {
        std::cerr  << "Could not create »" << filename << "«" << std::endl;
        return;
    }
    
    
    float* m_vertices = V;
    float* m_vertexNormals = Result_Normals;
    
    
    bool vertex_normal     = false;
    
    /* Add vertex element. */
    if ( m_vertices )
    {
        ply_add_element( oply, "vertex", m_numVertices );

        /* Add vertex properties: x, y, z, (r, g, b) */
        ply_add_scalar_property( oply, "x", PLY_FLOAT );
        ply_add_scalar_property( oply, "y", PLY_FLOAT );
        ply_add_scalar_property( oply, "z", PLY_FLOAT );
        
        /* Add normals if there are any. */
        if ( m_vertexNormals )
        {
			
            if ( m_numVertexNormals != m_numVertices )
            {
                std::cout << "Amount of vertices and normals"
                    << " does not match. Normals won't be written." << std::endl;
            }
            else
            {
                ply_add_scalar_property( oply, "nx", PLY_FLOAT );
                ply_add_scalar_property( oply, "ny", PLY_FLOAT );
                ply_add_scalar_property( oply, "nz", PLY_FLOAT );
                vertex_normal = true;
            }
        }
    }
    
    /* Write header to file. */
    if ( !ply_write_header( oply ) )
    {
        std::cerr  << "Could not write header." << std::endl;
        return;
    }
    
    int no_normal_index = 0;
    
    for (size_t i = 0; i < m_numVertices; i++ )
    {
		
		
        ply_write( oply, (double) m_vertices[ i * 3     ] ); /* x */
        ply_write( oply, (double) m_vertices[ i * 3 + 1 ] ); /* y */
        ply_write( oply, (double) m_vertices[ i * 3 + 2 ] ); /* z */
        
        if ( vertex_normal )
        {
			if(m_vertexNormals[i*3] == 0.0 && m_vertexNormals[i*3+1] == 0.0 && m_vertexNormals[i*3+2] == 0.0){
				no_normal_index = i;
			}
			
            ply_write( oply, (double) m_vertexNormals[ i * 3     ] ); /* nx */
            ply_write( oply, (double) m_vertexNormals[ i * 3 + 1 ] ); /* ny */
            ply_write( oply, (double) m_vertexNormals[ i * 3 + 2 ] ); /* nz */
        }
    }
    
    //~ std::cout << "last No normal index: " << no_normal_index << std::endl;
    
    if ( !ply_close( oply ) )
    {
       std::cerr  << "Could not close file." << std::endl;
    }

}

int main(int argc, char** argv){
	cuda_normals::Options opt(argc, argv);
    cout << opt << endl;
    
    ModelPtr model = ModelFactory::readModel(opt.inputFile());
    
    
    PointArray points;
	PointArray normals;
	
	normals.dim = 3;
	points.dim = 3;
	
	size_t num_points;
	
	if (model->m_pointCloud ) {
		
		points.elements = model->m_pointCloud->getPointArray(num_points).get();
		
	} else if (model->m_mesh) {
		
		points.elements = model->m_mesh->getVertexArray(num_points).get();
		
	}
	
	
	points.width = (int)num_points;
	
	
	floatArr normalsPtr = floatArr(new float[ num_points * 3 ]);
	normals.elements = normalsPtr.get();
	
	
	
	std::cout << "Constructing kd-tree..." << std::endl;
	CalcNormalsCuda calculator(points);
    std::cout << "Finished kd-tree construction." << std::endl;
	
	calculator.setK(50);
	calculator.setFlippoint(100000.0, 100000.0, 100000.0);
	calculator.setMethod("PCA");
	
    std::cout << "Start Normal Calculation..." << std::endl;
	calculator.start();
    std::cout << "Finished Normal Calculation." << std::endl;
	
	calculator.getNormals(normals);
	
	
    ModelPtr out_model;
    
    if(model->m_pointCloud) {
		
		PointBufferPtr buffer = model->m_pointCloud;
		buffer->setPointNormalArray( normalsPtr, normals.width);
		out_model = ModelPtr(new Model(buffer));
		
	} else if(model->m_mesh) {
		
		MeshBufferPtr buffer = model->m_mesh;
		buffer->setVertexNormalArray( normalsPtr, normals.width);
		out_model = ModelPtr(new Model(buffer));
		
	}
    

    ModelFactory::saveModel(out_model, "normals.ply");
    
}

int main2(int argc, char** argv) {
	
	
	
    const char* in_file = "/home/amock/datasets/polizei/raw/polizei30M_cut.ply";
    const char* out_file = "output_mesh.ply";
	
	size_t point_size;
	int point_dim = 3;
	
	PointArray points;
	PointArray normals;
	normals.dim = 3;
	
	
	
    std::cout << "Reading file " << in_file << " ..." << std::endl;
	readPlyFile(points, point_size, in_file, 1.0);
    std::cout << "Finished reading file. " << point_size << std::endl;

	points.width = point_size;
	points.dim = 3;
	
    std::cout << "Constructing kd-tree..." << std::endl;
	CalcNormalsCuda calculator(points);
    std::cout << "Finished kd-tree construction." << std::endl;
	
	calculator.setK(50);
	calculator.setFlippoint(100000.0, 100000.0, 100000.0);
	calculator.setMethod("PCA");
	
    std::cout << "Start Normal Calculation..." << std::endl;
	calculator.start();
    std::cout << "Finished Normal Calculation." << std::endl;
	
	calculator.getNormals(normals);
	
	
	
    std::cout << "Writing normals.." << std::endl;
	writePlyFile(points.elements, point_size, normals.elements, point_size, out_file);
    std::cout << "Finished writing." << std::endl;

	free(points.elements);
	return 1;
}
