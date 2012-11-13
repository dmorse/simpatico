#ifndef UTIL_AVERAGE_CPP
#define UTIL_AVERAGE_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Average.h"         // class header
#include <util/format/Dbl.h>
#include <util/format/Int.h>

#include <math.h>

namespace Util
{

   /*
   * Default constructor.
   */
   Average::Average(int blockFactor)
    : AverageStage(blockFactor),
      ParamComposite(),
      descendants_(),
      blockSum_(0.0),
      iBlock_(0),
      nSamplePerBlock_(0)
   {
      setClassName("Average");
      // Register self as first "descendant" AverageStage.
      descendants_.push_back(this);
   }

   Average::~Average()
   {}

   /*
   * Reset all accumulators and counters to zero.
   */
   void Average::clear()
   {

      blockSum_ = 0;
      iBlock_   = 0;
      AverageStage::clear();
   }

   /*
   * Read nSamplePerBlock from file.
   */
   void Average::readParam(std::istream& in)
   {
      readBegin(in, "Average");
      read<int>(in, "nSamplePerBlock", nSamplePerBlock_);
      readEnd(in);
   }

   /*
   * Set nSamplePerBlock parameter.
   */
   void Average::setNSamplePerBlock(int nSamplePerBlock)
   {
      nSamplePerBlock_ = nSamplePerBlock;
   }

   /*
   * Add a sampled value to the ensemble.
   */
   void Average::sample(double value)
   {
      std::ostream* outFilePtr = 0;
      sample(value, outFilePtr);
   }

   /*
   * Add a sampled value to the ensemble (private method)
   *
   * If outFilePtr != 0, output block averages to outFilePtr
   */
   void Average::sample(double value, std::ostream* outFilePtr)
   {

      // Add value to parent class
      AverageStage::sample(value);

      // Process block average for output
      if (nSamplePerBlock_ && outFilePtr) {

         // Increment block accumulators
         blockSum_  += value;
         ++iBlock_;

         if (iBlock_ == nSamplePerBlock_) {

            // Output block average to file
            *outFilePtr << Dbl(blockSum_/double(nSamplePerBlock_))
                        << std::endl;

            // Reset block accumulators
            blockSum_ = 0.0;
            iBlock_  = 0;

         }
      }
   }

   /*
   * Output statistical properties to file
   */
   void Average::output(std::ostream& out)
   {

      // Find first stage (descending) with nSample >= 16
      AverageStage* ptr = 0;
      int n = descendants_.size();
      int i = n;
      int nSample = 1;
      while (nSample < 16 && i > 0) {
         --i;
         ptr = descendants_[i];
         nSample = ptr->nSample();
      }

      double error  = ptr->error();
      double sigma  = error/sqrt(2.0*double(nSample-1));
      double weight = 1.0/(sigma*sigma);
      double sum    = error*weight;
      double norm   = weight;
      double aveErr = error;
      double oldSig;

      // Find weighted average within plateau
      bool next = true;
      while (next && i > 0) {
         oldSig = sigma;
         --i;
         ptr = descendants_[i];
         error = ptr->error();
         if (fabs(error - aveErr) < 2.0*oldSig) {
            nSample = ptr->nSample();
            sigma  = error/sqrt(2.0*double(nSample-1));
            weight = 1.0/(sigma*sigma);
            sum   += error*weight;
            norm  += weight;
            aveErr = sum/norm;
         } else {
            next = false;
         }
      }

      out <<  "Average   " << Dbl(average())         
          <<  "  +- " << Dbl(aveErr, 9, 2) << std::endl;
      out <<  "Variance  " << Dbl(variance())        << std::endl;
      out <<  "Std Dev   " << Dbl(stdDeviation())    << std::endl;
      out <<  std::endl;

      out << "Hierarchichal Error Analysis:" << std::endl;
      int interval;
      for (int i = 0; i < n; ++i) {
         ptr = descendants_[i];
         error    = ptr->error();
         nSample  = ptr->nSample();
         interval = ptr->stageInterval();
         if (nSample >= 16) {
            out << Int(i) 
                << Int(interval) 
                << Dbl(error) 
                << Dbl(error/sqrt(double(nSample)))
                << Int(nSample) << std::endl;
         }
      }
      out << std::endl;

   }

   /*
   * Append pointer to a descendant to descendants_ array.
   */
   void Average::registerDescendant(AverageStage* descendantPtr)
   {  descendants_.push_back(descendantPtr); }

   #if 0

   /**
   * Backup statistical accumulators to file.
   */
   void Average::backup(FILE *file)
   {}

   /**
   * Restore statistical accumulators from file.
   */
   void Average::restore(FILE *file)
   {}

   #endif

}
#endif
