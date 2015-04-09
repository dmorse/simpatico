/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "BondTensorAutoCorr.h"
#include <ddMd/analyzers/AutoCorrAnalyzer.tpp>
#include <ddMd/storage/BondStorage.h>
#include <ddMd/storage/GroupIterator.h>
#include <util/boundary/Boundary.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   BondTensorAutoCorr::BondTensorAutoCorr(Simulation& simulation) 
    : AutoCorrAnalyzer<Tensor, double>(simulation),
      bondTensor_()
   {  setClassName("BondTensorAutoCorr"); }

   /*
   * Destructor.
   */
   BondTensorAutoCorr::~BondTensorAutoCorr() 
   { }

   /*
   * Compute the bond tensor (Call on all processors).
   */
   void BondTensorAutoCorr::computeData()
   {
      MPI::Intracomm& communicator = simulation().domain().communicator();
      BondStorage& storage = simulation().bondStorage();
      Boundary& boundary = simulation().boundary();

      Tensor localTensor;
      Vector dr;
      double rsq;
      GroupIterator<2> iter;
      Atom* atom0Ptr;
      Atom* atom1Ptr;
      int isLocal0, isLocal1, i, j;


      // Iterate over bonds
      localTensor.zero();
      storage.begin(iter);
      for ( ; iter.notEnd(); ++iter) {
         atom0Ptr = iter->atomPtr(0);
         atom1Ptr = iter->atomPtr(1);
         isLocal0 = !(atom0Ptr->isGhost());
         isLocal1 = !(atom1Ptr->isGhost());
         rsq = boundary.distanceSq(atom0Ptr->position(), 
                                   atom1Ptr->position(), dr);
         if (rsq > 1.0E-6) {
            dr /= sqrt(rsq);
            assert(isLocal0 || isLocal1);
            if (isLocal0 && isLocal1) {
               // If both atoms are local
               for (i = 0; i < Dimension; ++i) {
                  for (j = 0; j < Dimension; ++j) {
                     localTensor(i, j) += dr[i]*dr[j];
                  }
               }
            } else {
               // If one atoms is local and one is a ghost
               for (i = 0; i < Dimension; ++i) {
                  for (j = 0; j < Dimension; ++j) {
                     localTensor(i, j) += 0.5*dr[i]*dr[j];
                  }
               }
            }
         }
      }

      // Reduce partial sums from all processors, store on the master.
      bondTensor_.zero();
      communicator.Reduce(&localTensor(0,0), &bondTensor_(0,0), 
                          Dimension*Dimension, MPI::DOUBLE, MPI::SUM, 0);

      if (communicator.Get_rank() != 0) {
         bondTensor_.zero();
      }
   }

   /*
   * Return the current bond tensor value.
   */
   Tensor BondTensorAutoCorr::data() 
   {  
      // Remove trace
      double trace = 0.0;
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         trace += bondTensor_(i,i);
      }
      trace = trace/double(Dimension);
      for (i = 0; i < Dimension; ++i) {
         bondTensor_(i,i) -= trace;
      }

      // Scale traceless symmetric tensor   
      double factor = 1.0/sqrt(10.0*simulation().boundary().volume());
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            bondTensor_(i,j) *= factor;
         }
      }

      return bondTensor_;
   }  

}
