#ifndef UTIL_AUTOCORR_STAGE_TPP
#define UTIL_AUTOCORR_STAGE_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrStage.h"

#include <util/accumulators/setToZero.h>
#include <util/accumulators/product.h>
#include <util/format/Int.h>
#include <util/format/write.h>

#include <complex>
#include <math.h>

using std::complex;

namespace Util
{

   /*
   * Constructor for rootPtr AutoCorrStage, with stageId = 0.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::AutoCorrStage()
    : buffer_(),
      corr_(),
      nCorr_(),
      sum_(),
      bufferCapacity_(0),
      nSample_(0),
      blockSum_(),
      nBlockSample_(0),
      stageInterval_(1),
      childPtr_(0),
      rootPtr_(0),
      stageId_(0),
      maxStageId_(10),
      blockFactor_(1)
   {  
      rootPtr_ = this; 
      setToZero(sum_);
      setToZero(blockSum_);
   }

   /*
   * Private constructor for stages with stageId > 0.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::AutoCorrStage(long stageInterval, 
                                               int stageId, int maxStageId,
                                               AutoCorrStage<Data, Product>* rootPtr, 
                                               int blockFactor)
    : buffer_(),
      corr_(),
      nCorr_(),
      sum_(),
      bufferCapacity_(rootPtr->bufferCapacity()),
      nSample_(0),
      blockSum_(),
      nBlockSample_(0),
      stageInterval_(stageInterval),
      childPtr_(0),
      rootPtr_(rootPtr),
      stageId_(stageId),
      maxStageId_(maxStageId),
      blockFactor_(blockFactor)
   { 
      allocate(); 
      setToZero(sum_);
      setToZero(blockSum_);
   }

   /*
   * Destructor.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::~AutoCorrStage()
   {
      if (childPtr_) {
         delete childPtr_;
      }
   }

   /*
   * Set buffer capacity and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::setParam(int bufferCapacity, 
                                               int maxStageId, int blockFactor)
   {
      bufferCapacity_ = bufferCapacity;
      maxStageId_ = maxStageId;
      blockFactor_ = blockFactor;
      allocate();
   }

   /*
   * Set previously allocated to initial empty state.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::clear()
   {
      setToZero(sum_);
      setToZero(blockSum_);
      nSample_ = 0;
      nBlockSample_ = 0;
      if (bufferCapacity_ > 0) {
         for (int i=0; i < bufferCapacity_; ++i) {
            setToZero(corr_[i]);
            nCorr_[i] = 0;
         }
         buffer_.clear();
      }
   }

   /*
   * Add a sampled value to the ensemble.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::sample(Data value)
   {
      if (bufferCapacity_ <= 0) {
         bufferCapacity_ = 64;
         blockFactor_ = 2;
         maxStageId_ = 10;
         allocate();
      }

      // Increment global accumulators
      sum_ += value;
      buffer_.append(value);
      for (int i=0; i < buffer_.size(); ++i) {
         corr_[i] += product(buffer_[i], value);
         ++nCorr_[i];
      };
      ++nSample_;

      // Increment block accumulators
      blockSum_ += value;
      ++nBlockSample_;

      if (nBlockSample_ == blockFactor_) {
         if (stageId_ < maxStageId_) {
            if (!childPtr_) {
               long nextStageInterval = stageInterval_*blockFactor_;
               int  nextStageId = stageId_ + 1;
               childPtr_ = new AutoCorrStage(nextStageInterval, nextStageId, 
                                             maxStageId_, rootPtr_, blockFactor_);
               rootPtr_->registerDescendant(childPtr_);
            }
            // Add block average as first value in child sequence 
            blockSum_ /= double(blockFactor_);
            childPtr_->sample(blockSum_);
         }
         // Reset block accumulators
         setToZero(blockSum_);
         nBlockSample_ = 0;
      }

   }

   /*
   * Serialize this AutoCorr.
   */
   template <typename Data, typename Product>
   template <class Archive>
   void 
   AutoCorrStage<Data, Product>::serialize(Archive& ar,
                                           const unsigned int version)
   {
      ar & buffer_;
      ar & corr_;
      ar & nCorr_;
      ar & sum_;
      ar & bufferCapacity_;
      ar & nSample_;
      isValid();

      ar & blockSum_;
      ar & nBlockSample_;
      ar & blockFactor_;

      // Constructor always sets stageInterval_ and stageId_

      // Does this stage have a child?
      int hasChild;
      if (Archive::is_saving()) {
         hasChild = (childPtr_ == 0) ? 0 : 1;
      }
      ar & hasChild;

      // Serialize child (if any)
      if (hasChild) {
         if (Archive::is_loading()) {
            long nextStageInterval = stageInterval_*blockFactor_;
            int  nextStageId       = stageId_ + 1;
            childPtr_ = new AutoCorrStage(nextStageInterval, 
                                          nextStageId,
                                          rootPtr_, blockFactor_);
            rootPtr_->registerDescendant(childPtr_);
         }
         ar & (*childPtr_);
      } else {
         if (Archive::is_loading()) {
            childPtr_ = 0;
         }
      }

   }

   /*
   * Return capacity of history buffer.
   */
   template <typename Data, typename Product>
   int AutoCorrStage<Data, Product>::bufferCapacity() const
   {  return bufferCapacity_; }

   /*
   * Return the number of sampled values.
   */
   template <typename Data, typename Product>
   long AutoCorrStage<Data, Product>::nSample() const
   {  return nSample_; }

   /*
   * Return the number of measured values per sample at this stage.
   */
   template <typename Data, typename Product>
   long AutoCorrStage<Data, Product>::stageInterval() const
   {  return stageInterval_; }

   /*
   * Return average of sampled values.
   */
   template <typename Data, typename Product>
   Data AutoCorrStage<Data, Product>::average() const
   { 
      return sum_/double(nSample_); 
   }

   /*
   * Calculate and output autocorrelation function.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::output(std::ostream& outFile)
   {
      //Data  ave = sum_/double(nSample_);
      //Product aveSq = product(ave, ave);

      int min;
      if (stageId_ == 0) {
         min = 0;
      } else {
         min = bufferCapacity_ / blockFactor_;
      }

      Product autocorr;
      for (int i = min; i < buffer_.size(); ++i) {
         autocorr = corr_[i]/double(nCorr_[i]);
         //autocorr = autocorr - aveSq;
         outFile << Int(i*stageInterval_) << " ";
         write<Product>(outFile, autocorr);
         outFile << std::endl;
      }
      if (childPtr_) {
         childPtr_->output(outFile);
      }
   }

   /*
   *  Return correlation time in unit of sampling interval
   */
   template <typename Data, typename Product>
   double AutoCorrStage<Data, Product>::corrTime() const
   {
      Data ave = sum_/double(nSample_);
      Product aveSq = product(ave, ave);
      Product variance = corr_[0]/double(nCorr_[0]);
      variance = variance - aveSq;
      Product autocorr, sum;
      setToZero(sum);
      int  size = buffer_.size();
      for (int i = 1; i < size/2; ++i) {
         autocorr = corr_[i]/double(nCorr_[i]);
         autocorr = autocorr - aveSq;
         sum += autocorr;
      }
      sum = sum/variance;
      return sum;
   }

   /*
   * Return autocorrelation at a given lag time index
   *
   * \param t the lag time index
   */
   template <typename Data, typename Product>
   Product AutoCorrStage<Data, Product>::autoCorrelation(int t) const
   {
      assert(t < buffer_.size());
      Product autocorr = corr_[t]/double(nCorr_[t]);
      Data ave  = sum_/double(nSample_);
      Product aveSq = product(ave, ave);
      autocorr -= aveSq;
      return autocorr;
   }

   // Private member functions

   /*
   * Allocate arrays and CyclicBuffer, and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::allocate()
   {
      if (bufferCapacity_ > 0) {
         corr_.allocate(bufferCapacity_);
         nCorr_.allocate(bufferCapacity_);
         buffer_.allocate(bufferCapacity_);
      }
      AutoCorrStage<Data, Product>::clear();
   }

   /*
   * Are capacities consistent?
   */
   template <typename Data, typename Product>
   bool AutoCorrStage<Data, Product>::isValid()
   {
      bool valid = true;
      if (bufferCapacity_ != corr_.capacity()) valid = false;
      if (bufferCapacity_ != nCorr_.capacity()) valid = false;
      if (bufferCapacity_ != buffer_.capacity()) valid = false;
      if (!valid) {
         UTIL_THROW("Invalid AutoCorr");
      }
      return valid;
   }

   template <typename Data, typename Product>
   void 
   AutoCorrStage<Data, Product>::registerDescendant(AutoCorrStage<Data, Product>* ptr)
   {}

}
#endif
