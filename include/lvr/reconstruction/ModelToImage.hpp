/* Copyright (C) 2016 Uni Osnabrück
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

/*
 * ModelToImage.h
 *
 *  Created on: Jan 25, 2017
 *      Author: Thomas Wiemann (twiemann@uos.de)
 */

#ifndef SRC_LIBLVR_RECONSTRUCTION_MODELTOIMAGE_HPP_
#define SRC_LIBLVR_RECONSTRUCTION_MODELTOIMAGE_HPP_

#include <lvr/io/Model.hpp>
#include <lvr/reconstruction/Projection.hpp>

#include <opencv/cv.h>
#include <algorithm>
#include <vector>

using std::vector;
using std::tuple;

namespace lvr {


class ModelToImage {
public:

    /// Pixelcoordinates with depth value
    typedef struct
    {
        int i;
        int j;
        float depth;
    } DepthPixel;

    typedef struct
    {
        float x;
        float y;
        float z;
    } PanoramaPoint;

    /// Image with single depth information
    typedef vector<vector<float> > DepthImage;

    /// Image with list of projected points at each pixel
    typedef vector<vector<PanoramaPoint> > PointListImage;



    ///
    /// \brief The ProjectionType enum
    ///
    enum ProjectionType {
        CYLINDRICAL, CONICAL, EQUALAREACYLINDRICAL,
        RECTILINEAR, PANNINI, STEREOGRAPHIC,
        ZAXIS, AZIMUTHAL
    };


    ///
    /// \brief Constructor
    /// \param buffer               A point buffer containing a point cloud
    /// \param projection           Type of used projection
    /// \param width                Desired image with. May be changed to fit
    ///                             panorama angles
    /// \param height               Desired image height. May be changed to fit
    ///                             panorama angles
    /// \param minZ                 Minimum depth value
    /// \param maxZ                 Maximum depth value
    /// \param minHorizontenAngle   Start of horizontal field of view in degrees
    /// \param maxHorizontalAngle   End of horizontal field of view in degrees
    /// \param mainVerticalAngle    Start of vertical field of view in degrees
    /// \param maxVerticalAngle     End of vertical field of view in degrees
    /// \param imageOptimization    If true, the aspect ration will be adapted to
    ///                             the chosen field of view
    /// \param leftHandedInputData  Set this to true of the scan points are in a
    ///                             left handed coordinate system (like 3dtk)
    ///
	ModelToImage(
			PointBufferPtr buffer,
            ProjectionType projection,
			int width, int height,
            int minZ, int maxZ,
			int minHorizontenAngle, int maxHorizontalAngle,
			int mainVerticalAngle, int maxVerticalAngle,
            bool imageOptimization,
            bool leftHandedInputData);

    ///
    /// \brief Writes the scan panaroma to an pgm file
    ///
    /// \param filename     Filename of the bitmap
    /// \param cutoff       Max range cutoff. Reduce this to enhance contrast on
    ///                     pixels with low depths.
    ///
    void writePGM(string filename, float cutoff);

    /// Destructor
	virtual ~ModelToImage();

    //// Returns an OpenCV image representation of the panarama
	void getCVMatrix(cv::Mat& image);

private:


    /// Pointer to projection
    Projection*     m_projection;

	/// Pointer to the initial point cloud
	PointBufferPtr 	m_points;

	/// Image width
	int				m_width;

	/// Image height
	int				m_height;

	/// Min horizontal opening angle
	int				m_minHAngle;

	/// Max horizontal opening angle
	int				m_maxHAngle;

	/// Min horizontal opening angle
	int				m_minVAngle;

	/// Max horizontal opening angle
	int				m_maxVAngle;

	/// Image optimization flag
	bool			m_optimize;

    int             m_xSize;

    int             m_ySize;

    int             m_maxWidth;

    int             m_maxHeight;

    int             m_minHeight;

    float           m_xFactor;

    float           m_yFactor;

    bool            m_leftHanded;
};

} /* namespace lvr */

#endif /* SRC_LIBLVR_RECONSTRUCTION_MODELTOIMAGE_HPP_ */