/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "WangLandauOutput.h"

#include <mcMd/simulation/Simulation.h>
#include <mcMd/mcMoves/semigrand/WangLandauMove.h>
#include <util/format/Dbl.h>
#include <util/format/Int.h>
#include <util/misc/FileMaster.h>

namespace McMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   WangLandauOutput::WangLandauOutput(McSystem& system)
    : SystemAnalyzer<McSystem>(system),
      outputFile_()
   {  setClassName("WangLandauOutput"); }

   /*
   * Read parameters from file and initialize. 
   */
   void WangLandauOutput::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);
      fileMaster().openOutputFile(outputFileName(".dat"), outputFile_);
   }

   /*
   * Load state from an archive.
   */
   void WangLandauOutput::loadParameters(Serializable::IArchive& ar)
   {  
      Analyzer::loadParameters(ar);
   }

   /*
   * Save state to an archive.
   */
   void WangLandauOutput::save(Serializable::OArchive& ar)
   {  ar & *this; }


   /* 
   * Sample occupation.
   */
   void WangLandauOutput::sample(long iStep) 
   {
      if (isAtInterval(iStep)) {
   //      distribution_.sample(mutatorPtr_->stateOccupancy(0));
   //      outputFile_ << iStep << "     " << mutatorPtr_->stateOccupancy(0) << std::endl;
           DArray<double> weights = getWeights();
      }
   }
 
   /*
   * Open and write summary output file.
   */
   void WangLandauOutput::output() 
   {
      //Close *.dat file
      outputFile_.close();
      //Open and write .hist file
   }
   
}
