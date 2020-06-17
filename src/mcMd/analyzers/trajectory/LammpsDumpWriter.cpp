/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsDumpWriter.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   LammpsDumpWriter::LammpsDumpWriter(System& system) 
    : TrajectoryWriter(system)
   {  setClassName("LammpsDumpWriter"); }

   /*
   * Destructor.
   */
   LammpsDumpWriter::~LammpsDumpWriter()
   {}
   
   /*
   * Write data that should appear once, at beginning of the file.
   */
   void LammpsDumpWriter::writeHeader()
   {};

   /*
   * Write data that should appear in every frame.
   */
   void LammpsDumpWriter::writeFrame(long iStep)
   {  system().writeConfig(outputFile_); }

}
