/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "OutputPairEnergies.h"
#include <util/misc/FileMaster.h>
#include <util/misc/ioUtil.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <sstream>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   OutputPairEnergies::OutputPairEnergies(Simulation& simulation) 
    : Analyzer(simulation),
      nSample_(0),
      isInitialized_(false)
   {  setClassName("OutputPairEnergies"); }

   /*
   * Read interval and outputFileName. 
   */
   void OutputPairEnergies::readParameters(std::istream& in) 
   {
      readInterval(in);
      readOutputFileName(in);

      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      isInitialized_ = true;
   }

   /*
   * Load internal state from an archive.
   */
   void OutputPairEnergies::loadParameters(Serializable::IArchive &ar)
   {
      loadInterval(ar);
      loadOutputFileName(ar);

      std::string filename;
      filename  = outputFileName();
      simulation().fileMaster().openOutputFile(filename, outputFile_);
      isInitialized_ = true;
   }

   /*
   * Save internal state to an archive.
   */
   void OutputPairEnergies::save(Serializable::OArchive &ar)
   {
      saveInterval(ar);
      saveOutputFileName(ar);
   }

  

   /*
   * Reset nSample_ counter.
   */
   void OutputPairEnergies::clear() 
   {  nSample_ = 0; }

   /*
   * Compute and output pair energies at regular intervals.
   */
   void OutputPairEnergies::sample(long iStep) 
   {
      if (isAtInterval(iStep))  {
         Simulation& sys = simulation();
         sys.computePairEnergies();
         if (sys.domain().isMaster()) {
            DMatrix<double> pair = sys.pairEnergies();
            for (int i = 0; i < simulation().nAtomType(); ++i){
               for (int j = 0; j < simulation().nAtomType(); ++j){
                  pair(i,j) = 0.5*( pair(i,j)+pair(j,i) );
                  pair(j,i) = pair(i,j);
               }
            }
            outputFile_ << Int(iStep, 10);
            for (int i = 0; i < simulation().nAtomType(); ++i){
               for (int j = 0; j < simulation().nAtomType(); ++j){
                  outputFile_ << Dbl(pair(i,j), 20);
               }
            }
            outputFile_  << std::endl;
         }

         ++nSample_;
      }
   }

}
