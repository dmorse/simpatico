#ifndef SPAN_TRAJECTORY_READER_CPP
#define SPAN_TRAJECTORY_READER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TrajectoryReader.h"

namespace SpAn
{

   using namespace Util;

   /*
   * Constructor.
   */
   TrajectoryReader::TrajectoryReader()
    : configurationPtr_(0)
   {  setClassName("TrajectoryReader"); }

   /*
   * Constructor.
   */
   TrajectoryReader::TrajectoryReader(Configuration& configuration)
    : configurationPtr_(&configuration)
   {  setClassName("TrajectoryReader"); }

   /*
   * Destructor.
   */
   TrajectoryReader::~TrajectoryReader()
   {}

}
#endif
