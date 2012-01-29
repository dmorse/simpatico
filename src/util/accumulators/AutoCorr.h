#ifndef AUTO_CORR_H
#define AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>    // Base class
#include <util/containers/RingBuffer.h>   // member
#include <util/containers/DArray.h>       // member

// Needed for implementation
#include <util/accumulators/setToZero.h>
#include <util/accumulators/product.h>
#include <util/space/Vector.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>
#include <util/format/write.h>

#include <complex>

using std::complex;

namespace Util
{

   /**
   * Calculates an auto-correlation function for a sequence of Data values.
   *   
   * This class calculates an autocorrelation function for a sequence
   * x(i) of values of a variable or object of type Data. The resulting
   * autocorrelation function is and array of values of type Product,
   * where C(j) = <x(i-j), x(i)>. Here <A,B> denotes an inner product of
   * type Product for objects A and B of type Data.
   *
   * The meaning of the inner product is defined for various data types b
   * the overloaded function Product product(Data, Data) that is defined
   * for double, complex and Vector data in the product.h file.
   *
   * The zero value for variables of type Data is returned by the overloaded
   * function void setToZero(Data) method defined in the setToData.h file.
   *
   * \ingroup Accumulators_Module
   */
   template <typename Data, typename Product>
   class AutoCorr : public ParamComposite
   {
   
   public:
  
      /**
      * Default constructor
      */
      AutoCorr();
   
      /**
      * Reset to empty state.
      */
      void clear();
   
      /**
      * Read buffer capacity, allocate memory and initialize.
      *
      * \param in input parameter stream.
      */
      void readParam(std::istream& in);
   
      /**
      * Set buffer capacity, allocate memory and initialize.
      *
      * \param bufferCapacity maximum number of values in history buffer.
      */
      void setParam(int bufferCapacity);

      /**
      * Sample a value.
      *
      * \param value current value
      */
      void sample(Data value);
   
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
      * Pack array into a block of memory.
      *
      * \param current pointer to current position in write buffer
      * \param end     pointer to end (one char* past last) of buffer
      */
      void pack(char*& current, char* end) const;

      /**
      * Unpack array from a block of memory.
      *
      * \param current pointer to current position in read buffer
      * \param end     pointer to end (one char* past last) of buffer
      */
      void unpack(char*& current, char* end);
       
      /**
      * Return packed size of this AutoCorr
      *
      * \return required sizeof packed buffer for this AutoCorr, in bytes.
      */
      int packedSize() const;

      /**
      * Serial to or from an Archive.
      * 
      * \param ar      input or output archive 
      * \param version id for file version
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:
   
      // Ring buffer containing a sequence of stored Data values.
      RingBuffer<Data>  buffer_; 
   
      // Array in which corr[j] = sum of values of <x(i-j), x(i)>
      DArray<Product>  corr_;
   
      // Array in which nCorr[i] = number of values added to corr[i]
      DArray<int>      nCorr_;
   
      // Sum of all previous values of x(t)
      Data        sum_;
   
      // Physical capacity (# of elements) of buffer, corr, and nCorr
      int         bufferCapacity_;
      
      // Total number of previous values of x(t)
      int         nSample_;
   
      /**
      * Allocate memory and initialize to empty state.
      */
      void allocate();

      /**
      * Is the internal state valid?
      */
      bool isValid();

   };

   /*
   * Default constructor.
   */
   template <typename Data, typename Product>
   AutoCorr<Data, Product>::AutoCorr() 
    : buffer_(),
      corr_(),
      nCorr_(),
      bufferCapacity_(0),
      nSample_(0)
   { setToZero(sum_); }
   
   /*
   * Read buffer capacity and allocate all required memory.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::readParam(std::istream& in)
   {
      readBegin(in, "AutoCorr");

      read<int>(in, "capacity", bufferCapacity_);
      allocate();

      readEnd(in);
   }
   
   /*
   * Set buffer capacity and initialize.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::setParam(int bufferCapacity)
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

      if (bufferCapacity_ > 0) {
         for (int i=0; i < bufferCapacity_; ++i) {
            setToZero(corr_[i]);
            nCorr_[i] = 0;
         }
         buffer_.clear();
      }
   }
   
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
   
   /*
   * Sample a single value from a time sequence.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::sample(Data value)
   {
      ++nSample_;
      sum_ += value;
      buffer_.append(value);
      for (int i=0; i < buffer_.size(); ++i) {
         corr_[i] += product(buffer_[i], value);
         ++nCorr_[i];
      };
   }
   
   /*
   * Return capacity of history buffer.
   */
   template <typename Data, typename Product>
   int AutoCorr<Data, Product>::bufferCapacity() const
   { return bufferCapacity_; }

   /*
   * Return number of values sampled thus far.
   */
   template <typename Data, typename Product>
   int AutoCorr<Data, Product>::nSample() const
   { return nSample_; }

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
         autocorr = autocorr - aveSq;
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
  
   /**
   * Return autocorrelation at a given lag time
   * 
   * \param t the lag time
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

   /*
   * Pack this AutoCorr into a memory block.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::pack(char*& current, char* end) const
   {
      if (current + packedSize() > end) {
         UTIL_THROW("Attempted write past end of send buffer");
      }
      buffer_.pack(current, end);
      
      corr_.pack(current, end);
      
      nCorr_.pack(current, end);
      
      Data* ptr = (Data *)(current);
      *ptr = sum_;
      current = (char*)(ptr + 1);
      
      int* ptr0 = (int*)(current);
      *ptr0 = bufferCapacity_;
      ++ptr0;
      *ptr0 = nSample_;
      ++ptr0;
      current = (char*)(ptr0);
   }

   /*
   * Unpack this AutoCorr from a memory block.
   */
   template <typename Data, typename Product>
   void AutoCorr<Data, Product>::unpack(char*& current, char* end) 
   {
      if (current + packedSize() > end) {
         UTIL_THROW("Attempted read past end of recv buffer");
      }
      buffer_.unpack(current, end);
      
      corr_.unpack(current, end);
      
      nCorr_.unpack(current, end);
      
      Data* ptr = (Data *)(current);
      sum_ = *ptr;
      current = (char*)(ptr + 1);
      
      int* ptr0 = (int*)(current);
      bufferCapacity_ = *ptr0;
      ++ptr0;
      nSample_ = *ptr0;
      ++ptr0;
      current = (char*)(ptr0);
   }

   /*
   * Return packed sizeof this AutoCorr, in bytes.
   */
   template <typename Data, typename Product>
   int AutoCorr<Data, Product>::packedSize() const
   {
      int size;
      size = buffer_.packedSize() + corr_.packedSize() + nCorr_.packedSize(); 
      size = size + sizeof(Data) + 2*sizeof(int);
      return size;
   }

   /*
   * Serialize this AutoCorr.
   */
   template <typename Data, typename Product>
   template <class Archive>
   void AutoCorr<Data, Product>::serialize(Archive& ar, 
                                           const unsigned int version)
   {
      ar & buffer_;
      ar & corr_;
      ar & nCorr_;
      ar & sum_;
      ar & bufferCapacity_;
      ar & nSample_;
      isValid();
   }

}
#endif
