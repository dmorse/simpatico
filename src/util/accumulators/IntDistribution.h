#ifndef UTIL_INT_DISTRIBUTION_H
#define UTIL_INT_DISTRIBUTION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/containers/DArray.h>

namespace Util
{

   /**
   * A distribution (or histogram) of values for an int variable.
   * 
   * \ingroup Accumulators_Module
   */
   class IntDistribution : public ParamComposite 
   {
   
   public:
  
      /** 
      * Default constructor
      */
      IntDistribution();
   
      /** 
      * Copy constructor.
      *
      * \param other object to be copied
      */
      IntDistribution(const IntDistribution& other);
   
      /** 
      * Assignment operator.
      *
      * \param other object to be assigned
      */
      IntDistribution& operator = (const IntDistribution& other);
   
      /** 
      * Destructor
      */
      virtual ~IntDistribution();
   
      /**
      * Read parameters from file and initialize.
      *
      * Read values of min, max, and nBin from file.
      * Allocate histogram array and clear all accumulators.
      *
      * \param in input parameter file stream
      */
      void readParam(std::istream& in);
  
      /** 
      * Set parameters and initialize.
      *
      * \param min  lower bound of range
      * \param max  upper bound of range
      */
      void setParam(int min, int max);
   
      /**
      * Clear (i.e., zero) previously allocated histogram.
      */
      void clear();
   
      /**
      * Sample a value.
      *
      * \param value current value
      */
      void sample(int value);
   
      /**
      * Output the distribution to file. 
      *
      * \param out output stream
      */
      void output(std::ostream& out);
  
      /**
      * Return the index of the bin for a value.
      *
      * \param value sampled value
      */
      int binIndex(int value);
  
      /** 
      * Get minimum value in range of histogram.
      */
      int min() const;
  
      /** 
      * Get maximum value in range of histogram.
      */
      int max() const;
    
      /** 
      * Get the number of bins
      */
      int nBin() const;
        
      /**
      * Serialize this Distribution to/from an archive.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   protected:
   
      DArray<long> histogram_;  ///< Histogram array.
      int   min_;               ///< minimum value.
      int   max_;               ///< maximum value.
      int   nBin_;              ///< number of bins.
      int   nSample_;           ///< Number of sampled values in Histogram.
      int   nReject_;           ///< Number of sampled values that were out of range.
   
   };

   // inline method definitions 
  
   /*
   * Return the index of the bin for a value.
   */
   inline int IntDistribution::binIndex(int value) 
   { return (value - min_); }
  
   /*
   * Get minimum value in range of histogram.
   */
   inline int IntDistribution::min() const
   {  return min_; }
  
   /*
   * Get maximum value in range of histogram.
   */
   inline int IntDistribution::max() const
   {  return max_; }
 
   /*
   * Get the number of bins
   */
   inline int IntDistribution::nBin() const
   {  return nBin_; }
        
   /*
   * Serialize this Distribution.
   */
   template <class Archive>
   void IntDistribution::serialize(Archive& ar, const unsigned int version)
   {
      ar & histogram_; 
      ar & min_;        
      ar & max_;    
      ar & nBin_;     
      ar & nSample_;   
      ar & nReject_;    
   }

}
#endif
