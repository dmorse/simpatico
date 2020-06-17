/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "ConfigDumpWriter.h"

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   ConfigDumpWriter::ConfigDumpWriter(System& system) 
    : TrajectoryWriter(system)
   {  setClassName("ConfigDumpWriter"); }

   /*
   * Destructor.
   */
   ConfigDumpWriter::~ConfigDumpWriter()
   {}
   
   /*
   * Write data that should appear once, at beginning of the file.
   */
   void ConfigDumpWriter::writeHeader()
   {};

   /*
   * Write data that should appear in every frame.
   */
   void ConfigDumpWriter::writeFrame(long iStep)
   {  system().writeConfig(outputFile_); }

}
