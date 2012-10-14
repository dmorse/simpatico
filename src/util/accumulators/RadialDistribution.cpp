#ifndef UTIL_RADIAL_DISTRIBUTION_CPP
#define UTIL_RADIAL_DISTRIBUTION_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "RadialDistribution.h"
#include <util/format/Int.h>
#include <util/format/Dbl.h>

namespace Util
{

   /* 
   * Default constructor.
   */
   RadialDistribution::RadialDistribution() 
    : Distribution(),
      norm_(0.0),
      nSnapshot_(0),
      outputIntegral_(false)
   {  setClassName("RadialDistribution"); }
   
   /* 
   * Copy constructor.
   */
   RadialDistribution::RadialDistribution(const RadialDistribution& other) 
    : Distribution(other),
      norm_(other.norm_),
      nSnapshot_(other.nSnapshot_),
      outputIntegral_(other.outputIntegral_)
   {}
   
   /* 
   * Assignment operator.
   */
   RadialDistribution& 
   RadialDistribution::operator = (const RadialDistribution& other) 
   {
      // Call base class assignment operator
      Distribution::operator=(other);

      // Assign the additional members
      norm_        = other.norm_;
      nSnapshot_   = other.nSnapshot_;

      return *this;
   }
   
   /* 
   * Read parameters and initialize.
   */
   void RadialDistribution::readParam(std::istream& in)
   {
      readBegin(in,"RadialDistribution");
      min_ = 0.0;
      read<double>(in, "max",  max_);
      read<int>(in, "nBin", nBin_);
      binWidth_ = (max_-min_)/double(nBin_);
      histogram_.allocate(nBin_);
      clear();
      readEnd(in);
   }

   /*
   * Set parameters and initialize.
   */
   void RadialDistribution::setParam(double max, int nBin)
   {
      min_   = 0.0;
      max_   = max;
      nBin_  = nBin;
      binWidth_  = (max_ - min_)/double(nBin_);
      histogram_.allocate(nBin_);
      clear();
   }
   
   /* 
   * Zero accumulators (virtual).
   */
   void RadialDistribution::clear()
   {
      Distribution::clear();
      nSnapshot_ = 0;
   }
   
   /* 
   * Set the factor used to normalize the RDF before output.
   */
   void RadialDistribution::setNorm(double norm)
   {  norm_ = norm; }
   
   /* 
   * Mark the beginning of a snapshot.
   */
   void RadialDistribution::beginSnapshot() 
   { ++nSnapshot_; }
   
   /* 
   * Set outputIntegral true/false to enable/disable output of spatial integral.
   */
   void RadialDistribution::setOutputIntegral(bool outputIntegral)
   { outputIntegral_ = outputIntegral; }
   
   /*
   * Output final results to file.
   */
   void RadialDistribution::output(std::ostream& out)
   {
      double r, rho, prefactor, dV, hist, integral;
      prefactor = 4.0*3.14159265359/3.0;
      prefactor = prefactor*binWidth_*binWidth_*binWidth_;
      integral  = 0.0;
      for (int i=0; i < nBin_; ++i) {
         r    = binWidth_*(double(i) + 0.5);
         dV   = prefactor*double(3*i*i + 3*i + 1);
         hist = double(histogram_[i])/double(nSnapshot_);
         rho  = hist/(dV*norm_);
         out << Dbl(r, 18, 8) << Dbl(rho, 18, 8);
	 if (outputIntegral_) {
            integral += (hist - dV*norm_)/norm_;
            out << Dbl(integral, 18, 8);
         }
         out <<  std::endl;
      }
   }
   
   #if 0
   /*
   * Backup statistical accumulators to file.
   */
   void RadialDistribution::backup(FILE *file) 
   {
      fprintf(file, "nSample       %i \n", nSample);
      fprintf(file, "nReject       %i \n", nReject);
      fprintf(file, "nSnapshot     %li \n", nSnapshot);
      for (int i=0; i < nBin; ++i) {
         fprintf(file, "%li ", histogram[i]);
      }
      fprintf(file, "\n");
   }
   
   /*
   * Restore statistical accumulators from file.
   */
   void RadialDistribution::restore(FILE *file) 
   {
      fscanf(file, "nSample       %i \n", &nSample);
      fscanf(file, "nReject       %i \n", &nReject);
      fscanf(file, "nSnapshot     %li \n", &nSnapshot);
      for (int i=0; i < nBin; ++i) {
         fscanf(file, "%li ", &histogram[i]);
      }
      fscanf(file, "\n");
   }
   #endif

}
#endif
