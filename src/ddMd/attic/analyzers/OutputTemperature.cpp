/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputTemperature.h"
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
   OutputTemperature::OutputTemperature(Simulation& simulation)
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("OutputTemperature"); }

   /*
   * Read interval and outputFileName.
   */
   void OutputTemperature::readParameters(std::istream& in)
   {
      readInterval(in);
      readOutputFileName(in);

      #if 0
      // Open output file
      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      #endif
      isInitialized_ = true;
   }


   /*
   * Load internal state from an archive.
   */
   void OutputTemperature::loadParameters(Serializable::IArchive &ar)
   {
      // Load parameter file parameters
      loadInterval(ar);
      loadOutputFileName(ar);

      // Load other data
      MpiLoader<Serializable::IArchive> loader(*this, ar);
      loader.load(nSample_);

      #if 0
      // Open output file
      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      #endif
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void OutputTemperature::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
      ar << nSample_;
   }

   /*
   * Clear nSample counter.
   */
   void OutputTemperature::clear()
   {  nSample_ = 0; }

   /*
   * Open outputfile
   */ 
   void OutputTemperature::setup()
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
   void OutputTemperature::sample(long iStep)
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computeKineticEnergy();
         simulation().atomStorage().computeNAtomTotal(simulation().domain().communicator());

         if (sys.domain().isMaster()) {
            double ndof = simulation().atomStorage().nAtomTotal()*3;
            double T_kinetic = sys.kineticEnergy()*2.0/ndof;
            outputFile_ << Int(iStep, 10)
                        << Dbl(T_kinetic, 20)
                        << std::endl;
         }

         ++nSample_;
      }
   }

}
