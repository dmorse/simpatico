/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "TensorAverage.h"         // class header
#include "Average.h"         
#include <util/format/Dbl.h>
#include <util/format/Int.h>

#include <math.h>

namespace Util
{

   /*
   * Default constructor.
   */
   TensorAverage::TensorAverage(int blockFactor)
    : ParamComposite(),
      iBlock_(0),
      nSamplePerBlock_(0)
   {
      setClassName("TensorAverage");
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            accumulators_[k].setBlockFactor(blockFactor);
            ++k;
         }
      }
   }

   /*
   * Destructor.
   */
   TensorAverage::~Average()
   {}

   /*
   * Reset all accumulators and counters to zero.
   */
   void TensorAverage::clear()
   {
      blockSum_ = 0;
      iBlock_   = 0;
      AverageStage::clear();
   }

   /*
   * Read nSamplePerBlock from file.
   */
   void TensorAverage::readParameters(std::istream& in)
   {  
      read<int>(in, "nSamplePerBlock", nSamplePerBlock_); 
      if (nSamplePerBlock_ < 0) {
         UTIL_THROW("Invalid input: nSamplePerBlock < 0");
      }
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            accumulators_[k].setNSamplePerBlock(nSamplePerBlock);
            ++k;
         }
      }
   }

   /*
   * Set nSamplePerBlock parameter.
   */
   void TensorAverage::setNSamplePerBlock(int nSamplePerBlock)
   {  
      if (nSamplePerBlock < 0) {
         UTIL_THROW("Attempt to set nSamplePerBlock < 0");
      }
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            accumulators_[k].setNSamplePerBlock(nSamplePerBlock);
            ++k;
         }
      }
   }

   /*
   * Load internal state from archive.
   */
   void TensorAverage::loadParameters(Serializable::IArchive &ar)
   {
      loadParameter<int>(ar, "nSamplePerBlock", nSamplePerBlock_); 
      if (nSamplePerBlock_ < 0) {
         UTIL_THROW("Loading value nSamplePerBlock < 0");
      }
      ar & iBlock_;
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            ar & accumulators_[k];
            ++k;
         }
      }
   }

   /*
   * Save internal state to archive.
   */
   void TensorAverage::save(Serializable::OArchive &ar)
   { 
      ar & nSamplePerBlock_;
      ar & iBlock_;
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            ar & accumulators_[k];
            ++k;
         }
      }
   }
   
   /*
   * Add a sampled value to the ensemble.
   */
   void TensorAverage::sample(const Tensor& value)
   {
      int i, j, k;
      k = 0;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            accumulators_[k].sample(value(i,j));
            ++k;
         }
      }
   }

   /*
   * Access accumulator associated with one component.
   */
   const Average& TensorAverage::accumulator(int i, int j)
   {
      int k = i*Dimension + j; 
      return accumulators_[k];
   }

}
