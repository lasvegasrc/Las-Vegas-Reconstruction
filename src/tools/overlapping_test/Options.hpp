/**
 * Copyright (C) 2013 Universität Osnabrück
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


#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <iostream>
#include <string>
#include <vector>
#include <boost/program_options.hpp>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::ostream;


namespace overlapping_test
{

using namespace boost::program_options;

/**
 * @brief A class to parse the program options for the reconstruction
 * executable.
 */
class Options
{
public:

	/**
	 * @brief Ctor. Parses the command parameters given to the main
	 *     function of the program
	 */
	Options(int argc, char** argv);
	virtual ~Options();

    string  outputFile() const
    {
        return (m_variables["outputFile"].as<string>());
    }

    string  inputFile() const
	{
        return (m_variables["inputFile"].as< vector<string> >())[0];
	}
	
	unsigned int maxLeafSize() const
	{
		return (m_variables["leafSize"].as<unsigned int>()); 
	}

	float overlapSize() const
	{
		return (m_variables["overlap"].as<float>());
	}

    bool blobMode() const
    {
        return (m_variables["blob"].as<bool>());
    }

    unsigned int nBlobs() const
    {
        return (m_variables["nblobs"].as<unsigned int>());
    }

private:

	/// The internally used variable map
	variables_map m_variables;

	/// The internally used option description
	options_description m_descr;

	/// The internally used positional option desription
	positional_options_description m_pdescr;

    float         m_overlap;
	unsigned int  m_leafsize;
    unsigned int  m_nBlobs;
    bool          m_blobMode;
    string        m_outputFile;
};

inline ostream& operator<<(ostream& os, const Options& o)
{
    os << "##### Settings: Overlap Tree Test #####" << endl;
	os << "Overlap Size: "<< o.overlapSize() << endl;
	os << "Max Leaf Size: " << o.maxLeafSize() << endl;
    
    return os;
}

} // namespace normals

#endif

