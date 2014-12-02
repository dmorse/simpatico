/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputBoxdim.h"
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/mpi/MpiLoader.h>
#include <util/misc/ioUtil.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   OutputBoxdim::OutputBoxdim(Simulation& simulation)
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("OutputBoxdim"); }

   /*
   * Read interval and outputFileName.
   */
   void OutputBoxdim::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);

      #if 0
      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      isInitialized_ = true;
      #endif
   }

   /*
   * Load internal state from an archive.
   */
   void OutputBoxdim::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void OutputBoxdim::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
   }

   /*
   * Reset nSample_
   */
   void OutputBoxdim::clear()
   {  nSample_ = 0; }

   /*
   * Open outputfile
   */ 
   void OutputBoxdim::setup()
   {
      if (simulation().domain().isMaster()) {
         std::string filename;
         filename  = outputFileName();
         simulation().fileMaster().openOutputFile(filename, outputFile_);
      }
   }

   /*
   * Dump configuration to file
   */
   void OutputBoxdim::sample(long iStep)
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeVirialStress();
         sys.computeKineticStress();
         if (sys.domain().isMaster()) {
            Vector L = sys.boundary().lengths();
            double V = sys.boundary().volume();
            outputFile_ << Int(iStep, 10)
                        << Dbl(L[0], 20)
                        << Dbl(L[1], 20)
                        << Dbl(L[2], 20)
                        << Dbl(V, 20)
                        << std::endl;
         }

         ++nSample_;
      }
   }

}
