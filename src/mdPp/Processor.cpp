#ifndef MDPP_PROCESSOR_CPP
#define MDPP_PROCESSOR_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Processor.h"

namespace MdPp 
{

    /**
    * Read parameters from file.
    */
    void Processor::readParameters(std::istream& in)
    {
       read<int>(in, "atomCapacity", atomCapacity_); 
    }

   
}
#endif
