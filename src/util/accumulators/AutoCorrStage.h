#ifndef UTIL_AUTOCORR_STAGE_H
#define UTIL_AUTOCORR_STAGE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // Base class
#include <util/containers/RingBuffer.h>   // member
#include <util/containers/DArray.h>       // member
#include <util/global.h>

namespace Util
{

   /**
   * Hierarchical auto-correlation function algorithm.
   *
   * This class calculates an autocorrelation function for a sequence
   * x(i) of values of a variable or object of type Data. The resulting
   * autocorrelation function is and array of values of type Product,
   * where C(j) = <x(i-j), x(i)>. Here <A,B> denotes an inner product 
   * of type Product for objects A and B of type Data.
   *
   * The meaning of the inner product is defined for various data 
   * types b the overloaded function Product product(Data, Data) that 
   * is defined for double, complex and Vector data in the product.h 
   * file.
   *
   * The zero value for variables of type Data is returned by the 
   * overloaded function void setToZero(Data) method defined in the 
   * setToData.h file.
   *
   * This class implements a hierarchical algorithm to calculate C(j).
   * The algorithm is implemented by a linked list of AutoCorrStage 
   * objects, in which each object in the list (except the first) 
   * calculates an autocorrelation for a sequence in which each value
   * in the sequence is the average of block of consecutive values in 
   * the primary sequence of measured values. Each object in this list
   * is assigned an integer chainId.  The "primary" AutoCorrStage 
   * object in this list, with chainId=0, calculates the autocorrelation 
   * for a "primary" sequence of values that are passed to the sample 
   * method of this object. For all n > 0, the object with chainId = n 
   * calculates the autocorrelation function for a sequence in which 
   * each value is an average of blockFactor values of the sequence 
   * managed by the parent object with chainId = n-1 or (equivalently) 
   * an average of blockFactor**n values of the primary sequence. New 
   * stages are added to this list dynamically as needed.
   *
   * \ingroup Accumulators_Module
   */
   template <typename Data, typename Product>
   class AutoCorrStage : public ParamComposite
   {

   public:

      /**
      * Constructor
      *
      * This constructor creates a primary AutoCorrStage object with
      * stageId = 0 and stageInterval = 1. A private constructor is
      * used to recursively create children of this object.
      *
      * \param blockFactor ratio of block sizes of subsequent stages
      */
      AutoCorrStage(int blockFactor = 4);

      /**
      * Destructor.
      *
      * Recursively destroy all descendant stages.
      */
      virtual ~AutoCorrStage();

      /**
      * Set the buffer (history) capacity.
      *
      * \param bufferCapacity max. number of values stored in buffer
      */
      void setCapacity(int bufferCapacity);

      /**
      * Sample a value.
      *
      * \param value current value
      */
      virtual void sample(Data value);

      /**
      * Register the creation of a descendant stage.
      *
      * This should be called only by a root stage.
      *
      * \param ptr pointer to a descendant AutoCorrStage.
      */
      virtual void registerDescendant(AutoCorrStage* ptr);

      /**
      * Serialize to/from an archive.
      *
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      //@}
      ///\name Accessors
      //@{

      /**
      * Output the autocorrelation function
      *
      * \param out output stream.
      */
      void output(std::ostream& out);

      /**
      * Return capacity of history buffer.
      */
      int bufferCapacity() const;

      /**
      * Return number of values sampled thus far.
      */
      int nSample() const;

      /**
      * Return the number of sampled values.
      */
      long nSample() const;

      /**
      * Return average of all sampled values.
      */
      Data average() const;

      /**
      * Numerical integration of autocorrelation function
      */
      double corrTime() const;

      /**
      * Return autocorrelation at a given lag time
      *
      * \param t the lag time
      */
      Product autoCorrelation(int t) const;

      /**
      * Return the number of sampled values per block at this stage.
      */
      long stageInterval() const;

      //@}

   protected:

      /**
      * Does this have a child AutoCorrStage?
      */
      bool hasChild() const;

      /**
      * Return the child AutoCorrStage by reference.
      */
      AutoCorrStage& child();

   private:

      // Ring buffer containing a sequence of stored Data values.
      RingBuffer<Data> buffer_;

      // Array in which corr[j] = sum of values of <x(i-j), x(i)>
      DArray<Product> corr_;

      // Array in which nCorr[i] = number of values added to corr[i]
      DArray<int> nCorr_;

      // Sum of all previous values of x(t)
      Data sum_;

      // Physical capacity (# of elements) of buffer, corr, and nCorr
      int bufferCapacity_;

      /// Number of sampled values.
      long nSample_;

      /// Sum of sampled values in the current block.
      double blockSum_;

      /// Number of values in the current block.
      long nBlockSample_;

      /// Number of measured values per sampled value at this stage.
      long stageInterval_;

      /// Pointer to child stage, if any.
      AutoCorrStage* childPtr_;

      /// Pointer to root stage. Null if this is the root stage.
      AutoCorrStage* rootPtr_;

      /// Stage index
      int stageId_;

      /// Number of samples per block.
      int blockFactor_;

      /**
      * Constructor for child objects.
      *
      * \param stageInterval  number of primary values per sample
      * \param stageId  integer id for this stage
      * \param rootPtr  pointer to root AutoCorrStage
      * \param blockFactor  ratio of block sizes of subsequent stages
      */
      AutoCorrStage(long stageInterval, int stageId, 
                    AutoCorrStage* rootPtr, int blockFactor);

      /**
      * Copy constructor - private and not implemented.
      */
      AutoCorrStage(const AutoCorrStage& other);

      /**
      * Assignment - private and not implemented.
      */
      AutoCorrStage& operator = (const AutoCorrStage& other);

      /**
      * Allocate memory and initialize to empty state.
      */
      void allocate();

      /**
      * Is the internal state valid?
      */
      bool isValid();

   };

   // Inline methods

   /*
   * Does this object have a child?  (protected)
   */
   inline bool AutoCorrStage::hasChild() const
   { return bool(childPtr_); }

   /*
   * Return child object by reference. (protected)
   */
   inline AutoCorrStage& AutoCorrStage::child()
   { return *childPtr_; }

}
#endif
