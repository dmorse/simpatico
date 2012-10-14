#ifndef UTIL_DISTRIBUTION_CPP
#define UTIL_DISTRIBUTION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Distribution.h"
#include <util/format/Dbl.h>
#include <util/global.h>

namespace Util
{

   /* 
   * Default constructor.
   */
   Distribution::Distribution() 
    : histogram_(),
      min_(0.0),
      max_(0.0),
      binWidth_(0.0),
      nBin_(0),
      nSample_(0),
      nReject_(0)
   {  setClassName("Distribution"); }
   
   /* 
   * Copy constructor.
   */
   Distribution::Distribution(const Distribution& other) 
    : histogram_(),
      min_(other.min_),
      max_(other.max_),
      binWidth_(other.binWidth_),
      nBin_(other.nBin_),
      nSample_(other.nSample_),
      nReject_(other.nReject_)
   {
      if (nBin_ > 0) {
         assert(other.histogram_.capacity() != 0);
         histogram_.allocate(nBin_);
         for (int i=0; i < nBin_; ++i) {
            histogram_[i] = other.histogram_[i];
         }
      } else {
         assert(nBin_ == 0);
         assert(histogram_.capacity() == 0);
         assert(nSample_ == 0);
         assert(nReject_ == 0);
      }
   }
   
   /* 
   * Assignment operator.
   */
   Distribution& Distribution::operator = (const Distribution& other) 
   {
      // Check for self assignment
      if (this == &other) return *this;

      // Check validity of other object
      if (other.nBin_ > 0) {
         assert(other.histogram_.capacity() != 0);
      } else {
         assert(other.nBin_ == 0);
         assert(other.histogram_.capacity() == 0);
         assert(other.nSample_ == 0);
         assert(other.nReject_ == 0);
      }

      // Assign primitive values
      min_      = other.min_;
      max_      = other.max_;
      binWidth_ = other.binWidth_;
      nBin_     = other.nBin_;
      nSample_  = other.nSample_;
      nReject_  = other.nReject_;

      // Allocate and copy histogram, if necessary
      if (nBin_ > 0) {
         histogram_.allocate(nBin_);
         for (int i=0; i < nBin_; ++i) {
            histogram_[i] = other.histogram_[i];
         }
      }

      return *this;
   }
   
   /* 
   * Destructor.
   */
   Distribution::~Distribution() 
   {}

   /* 
   * Read parameters and initialize.
   */
   void Distribution::readParam(std::istream& in)
   {
      readBegin(in,"Distribution");
      read<double>(in, "min", min_);
      read<double>(in, "max", max_);
      read<int>(in,   "nBin", nBin_);
      binWidth_  = (max_ - min_)/double(nBin_);
      histogram_.allocate(nBin_);
      clear();
      readEnd(in);
   }
  
   /*
   * Set parameters and initialize.
   *
   * \param min  lower bound of range
   * \param max  upper bound of range
   * \param nBin number of bins in range [min, max]
   */
   void Distribution::setParam(double min, double max, int nBin)
   {
      min_   = min;
      max_   = max;
      nBin_  = nBin;
      binWidth_  = (max_ - min_)/double(nBin_);
      histogram_.allocate(nBin_);
      clear();
   }  
   
   /* 
   * Zero all accumulators.
   */
   void Distribution::clear()
   {  
      nSample_ = 0; 
      nReject_ = 0; 
      for (int i=0; i < nBin_; ++i) {
         histogram_[i] = 0;
      }
   }
   
   /* 
   * Add a value to the histogram
   */
   void Distribution::sample(double value)
   {
      int i;
      if (value > min_ && value < max_) {
         i = binIndex(value);
         histogram_[i] += 1;
         nSample_ += 1;
      } else {
         nReject_ += 1;
      }
   }
   
   /* 
   * Output histogram
   */
   void Distribution::output(std::ostream& out) 
   {
      double x, rho;
      for (int i=0; i < nBin_; ++i) {
         x   = min_ + binWidth_*(double(i) + 0.5);
         rho = double(histogram_[i])/double(nSample_);
         rho = rho/binWidth_;
         out << Dbl(x, 18, 8) << Dbl(rho, 18, 8) << std::endl;
      }
   }
  
   #if 0 
   /* 
   *
   */
   void Distribution::backup(FILE *file) 
   {
      fprintf(file, "nSample       %i \n", nSample_);
      fprintf(file, "nReject       %i \n", nReject_);
      for (int i=0; i < nBin_; ++i) {
         fprintf(file, "%li ", histogram_[i]);
      }
      fprintf(file, "\n");
   }
   
   /* 
   *
   */
   void Distribution::restore(FILE *file) {
      fscanf(file, "nSample       %i \n", &nSample_);
      fscanf(file, "nReject       %i \n", &nReject_);
      for (int i=0; i < nBin_; ++i) {
         fscanf(file, "%li ", &histogram_[i]);
      }
      fscanf(file, "\n");
   }
   #endif

}
#endif
