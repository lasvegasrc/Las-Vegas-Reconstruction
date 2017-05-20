#ifndef IOUTILS_HPP
#define IOUTILS_HPP

#include <lvr/io/Timestamp.hpp>
#include <lvr/io/ModelFactory.hpp>
#include <boost/filesystem.hpp>

#include <Eigen/Dense>

#include <fstream>

namespace lvr
{

/**
 * @brief Returns a Eigen 4x4 maxtrix representation of the transformation
 *        represented in the given frame file.
 */
Eigen::Matrix4d getTransformationFromFrames(boost::filesystem::path& frames);

/**
 * @brief Transforms an slam6d transformation matrix into an Eigen 4x4 matrix.
 */
Eigen::Matrix4d buildTransformation(double* alignxf);

} // namespace lvr


#endif // IOUTILS_HPP

