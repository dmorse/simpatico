#ifndef DDMD_DOMAIN_CPP
#define DDMD_DOMAIN_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "Domain.h"
#include <util/space/Dimension.h>

namespace DdMd
{

   using namespace Util;

   /*
   * Constructor.
   */
   Domain::Domain()
    : Util::ParamComposite(),
      sourceRanks_(),
      destRanks_(),
      shift_(),
      gridDimensions_(),
      gridCoordinates_(),
      gridRank_(-1),
      gridIsPeriodic_(),
      #if UTIL_MPI
      intracommPtr_(0),
      #endif
      boundaryPtr_(0),
      isInitialized_(false)
   {  setClassName("Domain"); }

   /*
   * Destructor.
   */
   Domain::~Domain()
   {}

   #if UTIL_MPI
   /*
   * Set the grid intracommunicator.
   */
   void Domain::setGridCommunicator(MPI::Intracomm& intraCommunicator)
   {
      intracommPtr_ = &intraCommunicator;
      gridRank_ = intracommPtr_->Get_rank();
   }
   #else
   void Domain::setRank(int rank)
   {  gridRank_ = rank; }
   #endif

   /*
   * Set the Boundary object.
   */
   void Domain::setBoundary(Boundary& boundary)
   {  boundaryPtr_ = &boundary; }

   /*
   * Read parameters and initialize.
   */
   void Domain::readParameters(std::istream& in)
   {

      #ifdef UTIL_MPI
      if (intracommPtr_ == 0) {
         UTIL_THROW("Intra-communicator not set before readParam");
      } 
      #endif
      if (gridRank_ < 0) {
         UTIL_THROW("grid rank not set before readParam");
      }
      if (!hasBoundary()) {
         UTIL_THROW("Boundary not set before readParam");
      }
      if (isInitialized_) {
         UTIL_THROW("Already initialized");
      }

      // Read processor grid dimensions and initialize
      read<IntVector>(in, "gridDimensions", gridDimensions_);
      initialize();
   }
   
   /*
   * Read parameters and initialize.
   */
   void Domain::loadParameters(Serializable::IArchive& ar)
   {
      #ifdef UTIL_MPI
      if (intracommPtr_ == 0) {
         UTIL_THROW("Intra-communicator not set before readParam");
      } 
      #endif
      if (gridRank_ < 0) {
         UTIL_THROW("grid rank not set before readParam");
      }
      if (!hasBoundary()) {
         UTIL_THROW("Boundary not set before readParam");
      }
      if (isInitialized_) {
         UTIL_THROW("Already initialized");
      }

      // Read processor grid dimensions and initialize
      loadParameter<IntVector>(ar, "gridDimensions", gridDimensions_);
      initialize();
   }

   /*
   * Save internal state to an archive.
   */
   void Domain::save(Serializable::OArchive &ar)
   {  ar << gridDimensions_; }
  
   /*
   * Initialize data - called by readParameters and loadParameters (private).
   */
   void Domain::initialize() 
   {
      // Validate gridDimensions
      int nproc = 1;
      for (int i = 0; i < Dimension; i++) {
          if (gridDimensions_[i] <= 0) {
             UTIL_THROW("Processor grid dimensions must be greater than 0");
          }
          nproc *= gridDimensions_[i];
      }
      int commSize = 1;
      #ifdef UTIL_MPI
      commSize = intracommPtr_->Get_size();
      #endif
      if (nproc != commSize) {
         UTIL_THROW("Grid dimensions inconsistent with communicator size");
      }

      // Set grid dimensions
      grid_.setDimensions(gridDimensions_);

      // Mark all directions as periodic.
      for (int i = 0; i < Dimension; i++) {
         gridIsPeriodic_[i] = true;
      }

      // Find grid coordinates for this processor
      gridCoordinates_ = grid_.position(gridRank_);

      IntVector sourceCoordinates;
      int       i, j, k, jp;

      // Loop over Cartesian direction for send
      for (i = 0; i < Dimension; ++i) {

         // Loop over increment sign (+1 or -1 grid coordinate i)
         for (j = 0; j < 2; ++j) {

            // Increment to source along coordinate i
            k  = (j == 0 ? 1 : -1);
              
            // Complementary index (other direction)
            jp = (j == 0 ? 1 :  0);

            //Initialize coordinates of the partner processor
            sourceCoordinates    = gridCoordinates_;
            sourceCoordinates[i] = gridCoordinates_[i] + k;

            // Check for transfer past boundary
            if (gridIsPeriodic_[i] == true) {
               shift_(i, j) = grid_.shift(sourceCoordinates[i], i);
               sourceRanks_(i, j) = grid_.rank(sourceCoordinates);
               destRanks_(i, jp)  = sourceRanks_(i, j);
            } else {
               shift_(i, j)  = 0;
               if (!grid_.isInGrid(sourceCoordinates[i], i)) {
                  destRanks_(i, j)    = -1;
                  sourceRanks_(i, jp) = -1;
               }
            }
      
         }
      }
      isInitialized_ = true;
   }

   /*
   * Return one of the boundaries of the domain owned by this processor.
   */
   double Domain::domainBound(int i, int j) const
   {
      // Preconditions
      assert(isInitialized_);
      assert(boundaryPtr_);
      assert(i >= 0);
      assert(i < Dimension);
      assert(j >= 0);
      assert(j < 2);

      double dL;
      if (UTIL_ORTHOGONAL) {
         dL = boundaryPtr_->length(i) / double(gridDimensions_[i]);
      } else {
         dL = 1.0 / double(gridDimensions_[i]);
      }
      return (gridCoordinates_[i] + j)*dL;
   }

   /*
   * Return rank of processor whose domain contains a position Vector.
   */
   int Domain::ownerRank(const Vector& position) const
   {
      // Preconditions
      assert(isInitialized_);
      assert(boundaryPtr_);

      double dL;
      IntVector r;
      for (int i = 0; i < Dimension; ++i) {
         if (UTIL_ORTHOGONAL) {
            dL = boundaryPtr_->length(i) / double(gridDimensions_[i]);
         } else {
            dL = 1.0 / double(gridDimensions_[i]);
         }
         r[i] = int(position[i] / dL);
         if (r[i] < 0 || r[i] >= gridDimensions_[i]) {
            Log::file() << "Cart i   = " << i << std::endl;
            Log::file() << "position = " << position[i] << std::endl;
            Log::file() << "dL       = " << dL << std::endl;
            Log::file() << "r        = " << r[i] << std::endl;
            Log::file() << "gridDim  = " << gridDimensions_[i] << std::endl;
            UTIL_THROW("Invalid grid coordinate");
         }
      }
      return grid_.rank(r);
   }

   /*
   * Does a Vector lie in the domain owned by this processor?
   */
   bool Domain::isInDomain(const Vector& position) const
   {
      // Preconditions
      assert(isInitialized_);
      assert(boundaryPtr_);

      double dL;
      bool isIn = true;
      for (int i = 0; i < Dimension; ++i) {  
         if (UTIL_ORTHOGONAL) {
            dL = boundaryPtr_->length(i) / double(gridDimensions_[i]);
         } else {
            dL = 1.0 / double(gridDimensions_[i]);
         }
         if (position[i] <   gridCoordinates_[i]*dL) {
            isIn = false;
         }
         if (position[i] >= (gridCoordinates_[i] + 1)*dL) {
            isIn = false;
         }
      }
      return isIn;
   }

}
#endif
