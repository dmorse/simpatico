#ifndef UTIL_AVERAGE_H
#define UTIL_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/accumulators/AverageStage.h>  // base class
#include <util/param/ParamComposite.h>       // base class
#include <util/global.h>

#include <vector>

namespace Util
{

   /**
   * Calculates the average and variance of a sampled property.
   *
   * Average calculates block and global averages of a sampled value and 
   * its square, from which it obtains a global average and variance for
   * a sequence. A heirarchical blocking algorithm is used to estimate
   * the error on the average. No error estimate is provided for the
   * variance.
   *
   * The sample method also optionally outputs the block averages to file. 
   * The parameter nSamplePerBlock is the number of samples per block 
   * average. Setting this to zero suppresses output of block averages.
   *
   * \ingroup Accumulators_Module
   */
   class Average : public AverageStage, public ParamComposite
   {
   
   public:

      /**   
      * Constructor
      *
      * \param blockFactor ratio of block sizes for subsequent stages.
      */
      Average(int blockFactor = 2);
   
      /**   
      * Destructor
      */
      virtual ~Average();
   
      /**
      * Initialize accumulators to zero initial state.
      */
      void clear();
   
      /**
      * Read parameters from file and initialize.
      *
      * Read nSamplePerBlock. See setNSamplePerBlock() for
      * discussion of value.
      *
      * \param in input stream 
      */
      void readParam(std::istream& in);
   
      /**
      * Set nSamplePerBlock
      *
      * If nSamplePerBlock > 0, output a block average to file every
      * nSamplePerBlock samples.
      *
      * If nSamplePerBlock = 0, do not output block averages.
      *
      * \param nSamplePerBlock number of samples per block average output
      */
      void setNSamplePerBlock(int nSamplePerBlock);
  
      /**
      * Add a sampled value to the ensemble, and output block averages.
      *
      * \param value sampled value
      * \param out   output stream to which to write block averages
      */
      void sample(double value, std::ostream& out);
 
      /**
      * Add a sampled value to the ensemble.
      *
      * If blockFile != 0 and a datafile has been set, this method
      * outputs block averages to the datafile.
      *
      * \param value  sampled value
      */
      void sample(double value);
 
      /**
      * Output final statistical properties to file.
      *
      * \param out output stream
      */
      void output(std::ostream& out);
  
      /**
      * Get nSamplePerBlock. 
      *
      * A zero value indicates that no block averages will be output.
      */
      int nSamplePerBlock();

      /**
      * Add pointer to a descendant to an array.
      */
      virtual void registerDescendant(AverageStage* descendantPtr);
         
      /**
      * Serialize this Average to or from an archive.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:
   
      /**
      * Add a sampled value to the ensemble, and output block averages.
      *
      * Outputs block averages to stream *outPtr iff outPtr != 0.
      *
      * \param value  sampled value
      * \param outPtr pointer to output stream
      */
      void sample(double value, std::ostream* outPtr);

      /// Array of pointers to descendants.
      std::vector<AverageStage*> descendants_;
 
      /// Sum of values in current output block
      double blockSum_;    

      /// Number of samples in current output block.
      int    iBlock_;     
      
      /// Number of sampled values per output block.
      int    nSamplePerBlock_;  

      /// Private and not implemented to prohibit copying.
      Average(const Average& other);

      /// Private and not implemented to prohibit assignment.
      Average& operator = (const Average& other);
     
   };
   
   // Inline method definitions
 
   /*
   * Add a sampled value to the ensemble, and output block averages.
   */
   inline void Average::sample(double value, std::ostream& out)
   {  sample(value, &out); }
 
   /*
   * Get nSamplePerBlock. 
   */
   inline int Average::nSamplePerBlock()
   {  return nSamplePerBlock_; }

   /*
   * Serialize this Average.
   */
   template <class Archive>
   void Average::serialize(Archive& ar, const unsigned int version)
   {
      AverageStage::serialize(ar, version);
      ar & blockSum_;    
      ar & iBlock_;     
      ar & nSamplePerBlock_;  
   }
 
}
#endif
