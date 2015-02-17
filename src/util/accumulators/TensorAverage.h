#ifndef UTIL_TENSOR_AVERAGE_H
#define UTIL_TENSOR_AVERAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class
#include <util/accumulators/Average.h>  // member 
#include <util/global.h>

#include <vector>

namespace Util
{

   /**
   * Calculates the averages of components of a Tensor.
   *
   * TensorAverage calculates block and global averages of components
   * of a sampled Tensor.  A hierarchical blocking algorithm is used 
   * to estimate the error on the average of each component. 
   *
   * The sample function of also optionally calculates block averages,
   * which can be useful for reducing how frequently values are logged
   * to a file. The parameter nSamplePerBlock is the number of samples 
   * per block average. This is initialized to zero. A zero value 
   * disables calculation of block averages. An overloaded method of 
   * the sample function that takes an std::ostream file as an argument 
   * outputs block averages to file as blocks are completed.
   *
   * \ingroup Accumulators_Module
   */
   class TensorAverage : public ParamComposite
   {

   public:

      /**
      * Constructor
      *
      * \param blockFactor ratio of block sizes for subsequent stages.
      */
      TensorAverage(int blockFactor = 2);

      /**
      * Destructor
      */
      virtual ~TensorAverage();

      /**
      * Read parameter nSamplePerBlock from file and initialize.
      *
      * See setNSamplePerBlock() for discussion of value.
      *
      * \param in input stream
      */
      void readParameters(std::istream& in);

      /**
      * Set nSamplePerBlock.
      *
      * If nSamplePerBlock > 0, the sample function will increment block
      * averages, and reset the average every nSamplePerBlock samples.
      *
      * If nSamplePerBlock == 0, block averaging is disabled. This is the
      * default (i.e., the initial value set in the constructor).
      *
      * \param nSamplePerBlock number of samples per block average output
      */
      void setNSamplePerBlock(int nSamplePerBlock);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Serialize this Average to or from an archive.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /**
      * Clear all accumulators, set to empty initial state.
      */
      void clear();

      /**
      * Add a sampled value to the ensemble.
      *
      * \param value sampled value
      */
      void sample(double value);

      /**
      * Output final statistical properties to file.
      *
      * \param out output stream
      */
      void output(std::ostream& out);

      /**
      * Get number of samples per block average.
      *
      * A zero value indicates that block averaging is disabled.
      */
      int nSamplePerBlock() const;

      /**
      * Get number of samples in current block average.
      *
      * Return 0 if block averaging disabled, if !nSamplePerBlock.
      */
      int iBlock() const;

      /**
      * Is the current block average complete?
      * 
      * \return (iBlock > 0) && (iBlock == nSamplePerBlock)
      */
      bool isBlockComplete() const;
   
   private:

      FArray<Average, Dimension*Dimension> accumulators_;

      /// Number of samples in current output block.
      int iBlock_;

      /// Number of sampled values per output block.
      int nSamplePerBlock_;

      /// Private and not implemented to prohibit copying.
      TensorAverage(const TensorAverage& other);

      /// Private and not implemented to prohibit assignment.
      TensorAverage& operator = (const TensorAverage& other);

   };

   // Inline method definitions

   /*
   * Get nSamplePerBlock, number of samples per block average.
   */
   inline int Average::nSamplePerBlock() const
   {  return nSamplePerBlock_; }

   /*
   * Get iBlock, number of samples in current block average.
   */
   inline int Average::iBlock() const
   {  return iBlock_; }

   /*
   * Is the current block average complete?
   */
   inline bool Average::isBlockComplete() const
   { 
      return (iBlock_ && (iBlock_ == nSamplePerBlock_));
   }

   /*
   * Get current block average.
   */
   inline Tensor Average::blockAverage() const
   {
      Tensor average;
      if (iBlock_ == 0) {
         UTIL_THROW("Attempt to get block average with no data");
      }
      return blockSum_/double(iBlock_);
   }

   /*
   * Serialize this Average.
   */
   template <class Archive>
   void Average::serialize(Archive& ar, const unsigned int version)
   {
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            ar << accumulators_[k];
            ++k;
         }
      }
      ar & iBlock_;
      ar & nSamplePerBlock_;
   }

}
#endif
