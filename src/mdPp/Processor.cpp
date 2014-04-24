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
       read<int>(in, "bondCapacity", bondCapacity_); 
       // etc. for angles dihedrals

       atoms_.allocate(atomCapacity_);
       bonds_.allocate(bondCapacity_);
       // etc. for angles dihedrals

       //readParamComposite(in, "AnalyzerManager", AnalyzerManager_);

       read<std::string>(in, "configIoName", configIoName_);
       configIoPtr_ = configIoFactory_.factory(configIoName_);

       read<std::string>(in, "configFileName", configFileName_);
    }
   
}
#endif
