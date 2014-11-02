#ifndef UTIL_AUTOCORR_STAGE_TPP
#define UTIL_AUTOCORR_STAGE_TPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrStage.h"

#include <util/accumulators/setToZero.h>
#include <util/accumulators/product.h>
#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/format/write.h>

#include <complex>
#include <math.h>

using std::complex;

namespace Util
{

   /*
   * Default constructor.
   */
   template <typename Data, typename Product>
   AutoCorr<Data, Product>::AutoCorr()
   {
      setClassName("AutoCorr");
      setToZero(sum_);
   }

   /*
   * Constructor for rootPtr AutoCorrStage, with stageId = 0.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::AutoCorrStage(int blockFactor)
    : buffer_(),
      corr_(),
      nCorr_(),
      bufferCapacity_(0),
      nSample_(0)
      blockSum_(0.0),
      nBlockSample_(0),
      stageInterval_(1),
      childPtr_(0),
      rootPtr_(0),
      stageId_(0),
      blockFactor_(blockFactor)
   {  rootPtr_ = this; }

   /*
   * Constructor for dynamically generated objects with stageId > 0.
   */
   template <typename Data, typename Product>
   AutoCorrStage<Data, Product>::AutoCorrStage(long stageInterval, 
                                               int stageId, 
                                               AutoCorrStage* rootPtr, 
                                               int blockFactor)
    : buffer_(),
      corr_(),
      nCorr_(),
      bufferCapacity_(rootPtr->bufferCapacity()),
      nSample_(0)
      blockSum_(0.0),
      nBlockSample_(0),
      stageInterval_(stageInterval),
      childPtr_(0),
      rootPtr_(rootPtr),
      stageId_(stageId),
      blockFactor_(blockFactor)
   { 
      allocate(); 
   }

   /*
   * Destructor
   */
   template <typename Data, typename Product>
   AutoCorr<Data, Product>::~AutoCorr()
   {}

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
   void AutoCorr<Data, Product>::setCapacity(int bufferCapacity)
   {
      bufferCapacity_ = bufferCapacity;
      allocate();
   }

   /*
   * Set previously allocated to initial empty state.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::clear()
   {
      setToZero(sum_);
      nSample_ = 0;
      blockSum_ = 0.0;
      nBlockSample_ = 0;
      if (bufferCapacity_ > 0) {
         for (int i=0; i < bufferCapacity_; ++i) {
            setToZero(corr_[i]);
            nCorr_[i] = 0;
         }
         buffer_.clear();
      }
      if (childPtr_) {
         delete childPtr_;
      }
   }

   /*
   * Add a sampled value to the ensemble.
   */
   template <typename Data, typename Product>
   void AutoCorrStage<Data, Product>::sample(double value)
   {

      // Increment global accumulators
      ++nSample_;
      sum_ += value;
      buffer_.append(value);
      for (int i=0; i < buffer_.size(); ++i) {
         corr_[i] += product(buffer_[i], value);
         ++nCorr_[i];
      };

      // Increment block accumulators
      blockSum_ += value;
      ++nBlockSample_;

      if (nBlockSample_ == blockFactor_) {

         if (!childPtr_) {
            long nextStageInterval = stageInterval_*blockFactor_;
            int  nextStageId = stageId_ + 1;
            childPtr_ = new AutoCorrStage(nextStageInterval, nextStageId, 
                                         rootPtr_, blockFactor_);
            rootPtr_->registerDescendant(childPtr_);
         }

         // Add block average as first value in child sequence 
         childPtr_->sample(blockSum_ / double(blockFactor_));

         // Reset block accumulators
         blockSum_ = 0.0;
         nBlockSample_ = 0;

      }

   }

   /*
   * Serialize this AutoCorr.
   */
   template <typename Data, typename Product>
   template <class Archive>
   void AutoCorr<Data, Product>::serialize(Archive& ar,
                                           const unsigned int version)
   {
      ar & bufferCapacity_;
      ar & buffer_;
      ar & corr_;
      ar & nCorr_;
      ar & sum_;
      ar & nSample_;
      isValid();

      ar & blockSum_;
      ar & nBlockSample_;
      ar & stageInterval_;
      ar & blockFactor_;

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
   int AutoCorr<Data, Product>::bufferCapacity() const
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
   Data AutoCorr<Data, Product>::average() const
   {
      Data ave  = sum_;
      ave /= double(nSample_);
      return ave;
   }

   /*
   * Final output
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::output(std::ostream& outFile)
   {
      Data    ave;
      Product autocorr;
      Product aveSq;

      // Calculate and output average of sampled values
      ave  = sum_;
      ave /= double(nSample_);
      aveSq = product(ave, ave);

      // Calculate and output autocorrelation
      for (int i = 0; i < buffer_.size(); ++i) {
         autocorr = corr_[i]/double(nCorr_[i]);
         //autocorr = autocorr - aveSq;
         outFile << Int(i);
         write<Product>(outFile, autocorr);
         outFile << std::endl;
      }

   }

   /*
   *  Return correlation time in unit of sampling interval
   */
   template <typename Data, typename Product>
   double AutoCorr<Data, Product>::corrTime() const
   {
      Data    ave;
      Product aveSq, variance, autocorr, sum;
      int     size;

      size = buffer_.size();

      // Calculate average sampled values
      ave   = sum_/double(nSample_);
      ave  /= double(nSample_);
      aveSq = product(ave, ave);

      // Calculate variance sampled values
      variance = corr_[0]/double(nCorr_[0]);
      variance = variance - aveSq;

      setToZero(sum);
      for (int i = 1; i < size/2; ++i) {
         autocorr = corr_[i]/double(nCorr_[i]);
         autocorr = autocorr - aveSq;
         sum     += autocorr;
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
   Product AutoCorr<Data, Product>::autoCorrelation(int t) const
   {
      Data    ave;
      Product autocorr;
      Product aveSq;

      // Calculate average of sampled values
      ave  = sum_;
      ave /= double(nSample_);
      aveSq = product(ave, ave);

      // Calculate and return autocorrelation
      assert(t < buffer_.size());
      autocorr = corr_[t]/double(nCorr_[t]);
      autocorr = autocorr - aveSq;

      return autocorr;
   }

   // Private member functions

   /*
   * Allocate arrays and CyclicBuffer, and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::allocate()
   {
      if (bufferCapacity_ > 0) {
         corr_.allocate(bufferCapacity_);
         nCorr_.allocate(bufferCapacity_);
         buffer_.allocate(bufferCapacity_);
      }
      clear();
   }

   /*
   * Are capacities consistent?
   */
   template <typename Data, typename Product>
   bool AutoCorr<Data, Product>::isValid()
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
   AutoCorrStage<Data, Product>::registerDescendant(AutoCorrStage* ptr)
   {}

}
#endif
