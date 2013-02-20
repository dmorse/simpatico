#ifndef UTIL_AVERAGE_STAGE_H
#define UTIL_AVERAGE_STAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

namespace Util
{

   /**
   * A hierarchical block average algorithm.
   *
   * This class implements a hierarchical block averaging algorithm
   * that calculates the average and variance of a sequence of values,
   * and also calculates the variance for sequences of block averages,
   * for multiple block sizes. The algorithm is implemented by a
   * chain of AverageStage objects, in which each object is assigned
   * an integer chainId. The "primary" AverageStage object (chainId=0)
   * in this chain calculates the average and variance for a "primary"
   * sequence of measured values that are passed as parameters to
   * its sample method. This object hold a pointer to an AverageStage
   * with chainId=1 that calculates the variance of a sequence in which
   * each value is the average of a block of blockFactor consecutive
   * values in the primary sequence. The object with chainId=1 holds a
   * pointer to an object with chainId=2 that calculates the variance
   * of a sequence in which each value is the average of a block of
   * blockFactor values of the sequence processed by chainId=1, or
   * blockFactor**2 values of the primary sequence. In general, the
   * object with chainId=n calculates the variance of a sequence in
   * which each value is an average of blockFactor**n values of the
   * primary sequence.  The integer parameter blockFactor is passed
   * to the constructor of the primary AverageStage object as a
   * parameter, which is set to blockFactor=2 by default. New stages
   * are added to this chain of objects dynamically as the length of
   * the primary sequence grows: After the primary AverageStage has
   * sampled m values, there exists one AverageStage objects for each
   * value of 0 <= * stageId <= n, where n is the largest integer for
   * which blockFactor**n <= m.
   *
   * The outputError() method of the primary AverageStage outputs
   * a sequence of estimates of the error of the average that are
   * derived different stages of block averaging. The estimate
   * obtained from each stage is given by sqrt(variance/nSample),
   * where variance is the variance of the sequence of block averages
   * processed by that stage, and nSample is the number of such block
   * averages thus far. This estimate is accurate only when the block
   * averages are long enough to be uncorrelated. A reliable
   * estimate of the error of the average can be obtained if and
   * only if several stages of block averaging yield statistically
   * indistinguishable error estimates.
   *
   * Algorithm:
   *
   * The averaging algorithm is implemented by the sample() method.
   * Each time the sample method is called with a new value, the
   * method increments a sum of values and sum of squares,from which
   * the average and variance can be calculated. It also increments a
   * block sum that is reset to zero every blockFactor samples. Each
   * AverageStage object can hold a pointer to a child AverageStage
   * object. When an AverageStage is instantiated, that pointer is
   * null. After the first blockFactor samples, however, the sample
   * method creates a new child object. At this point, and every
   * blockFactor steps thereafter, the sample method of the parent
   * passes a block average to the sample method of its child and
   * then resets the block sum to zero. The resulting chain of objects
   * is extended as needed: After blockFactor*blockFactor samples have 
   * been passed to a parent, so that blockFactor block average values
   * have been passed to a child object, the child will create a
   * grandchild. The public interface of the primary AverageStage
   * object allows all of the objects in the resulting chain to
   * report error estimates, via the recursive outputError() method,
   * but does not allow any other form of access to descendants of
   * the primary AverageStage object.
   *
   * \ingroup Accumulators_Module
   */
   class AverageStage
   {

   public:

      /**
      * Constructor
      *
      * This constructor creates a primary AverageStage object with
      * stageId = 0 and stageInterval = 1. A private constructor is
      * used to recursively create children of this object.
      *
      * \param blockFactor ratio of block sizes of subsequent stages
      */
      AverageStage(int blockFactor = 2);

      /**
      * Destructor.
      *
      * Recursively destroy all children.
      */
      virtual ~AverageStage();

      /**
      * Initialize all accumulators and recursively destroy all children.
      */
      virtual void clear();

      /**
      * Register the creation of a descendant stage.
      *
      * This should be called only by a root stage.
      *
      * \param descendantPtr pointer to a descendant AverageStage.
      */
      virtual void registerDescendant(AverageStage* descendantPtr);

      /**
      * Add a sampled value to the ensemble.
      *
      * \param value sampled value
      */
      virtual void sample(double value);

      /**
      * Add a sampled value to the ensemble.
      *
      * \param ar       input or output archive
      * \param version  file version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      ///\name Accessors
      //@{

      /**
      * Return the average of all sampled values.
      */
      double average() const;

      /**
      * Return the variance of all sampled values.
      */
      double variance() const;

      /**
      * Return the standard deviation of all sampled values.
      *
      * \return sqrt(variance())
      */
      double stdDeviation() const;

      /**
      * Return an estimate for the std deviation of the average.
      *
      * \return sqrt(variance()/nSample())
      */
      double error() const;

      /**
      * Return the number of sampled values.
      */
      long nSample() const;

      /**
      * Return the number of sampled values per block at this stage.
      */
      long stageInterval() const;

      //@}

   protected:

      /**
      * Does this object have a child AverageStage for block averages?
      */
      bool hasChild() const;

      /**
      * Return the child AverageStage by reference.
      */
      AverageStage& child();

   private:

      /// Sum of all sampled values.
      double sum_;

      /// Sum of squares of all sampled values.
      double sumSq_;

      /// Sum of sampled values in the current block.
      double blockSum_;

      /// Number of sampled values.
      long   nSample_;

      /// Number of values in the current block.
      long   nBlockSample_;

      /// Number of measured values per sampled value at this stage.
      long   stageInterval_;

      /// Pointer to child stage, if any.
      AverageStage* childPtr_;

      /// Pointer to root stage. Null if this is the root stage.
      AverageStage* rootPtr_;

      /// Stage index
      int stageId_;

      /// Number of samples per block.
      int blockFactor_;

      /**
      * Constructor for child objects.
      *
      * \param stageInterval number of measured values per sample at this stage
      * \param stageId       integer id for this stage
      * \param rootPtr       pointer to root AverageStage
      * \param blockFactor   ratio of block sizes of subsequent stages
      */
      AverageStage(long stageInterval, int stageId, AverageStage* rootPtr, int blockFactor);

      /**
      * Copy constructor - private and not implemented.
      */
      AverageStage(const AverageStage& other);

      /**
      * Assignment - private and not implemented.
      */
      AverageStage& operator = (const AverageStage& other);

   };

   // Inline methods

   /*
   * Does this object have a child?  (protected)
   */
   inline bool AverageStage::hasChild() const
   { return bool(childPtr_); }

   /*
   * Return child object by reference. (protected)
   */
   inline AverageStage& AverageStage::child()
   { return *childPtr_; }

   // Method template

   /*
   * Serialize this stage.
   */
   template <class Archive>
   void AverageStage::serialize(Archive& ar, const unsigned int version)
   {
      ar & sum_;
      ar & sumSq_;
      ar & blockSum_;
      ar & nSample_;
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
            childPtr_ = new AverageStage(nextStageInterval, nextStageId,
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

}
#endif
