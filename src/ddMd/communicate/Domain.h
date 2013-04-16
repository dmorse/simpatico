#ifndef DDMD_DOMAIN_H
#define DDMD_DOMAIN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class
#include <util/boundary/Boundary.h>     // typedef used in interface
#include <util/containers/FMatrix.h>    // member template
#include <util/containers/FArray.h>     // member template
#include <util/space/IntVector.h>        // member
#include <util/space/Grid.h>             // member
#include <util/space/Dimension.h>        // constant expression
#include <util/global.h>

namespace DdMd
{

   using namespace Util;

   /**
   * Decomposition of the system into domains associated with processors.
   * 
   * \ingroup DdMd_Communicate_Module
   */
   class Domain : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Domain();

      /**
      * Destructor.
      */
      virtual ~Domain();

      // Mutators

      #if UTIL_MPI

      /** 
      * Set the grid communicator.
      *
      * \param communicator Intramcommunicator to be used.
      */
      void setGridCommunicator(MPI::Intracomm& communicator);

      #else

      void setRank(int Rank);

      #endif

      /** 
      * Set the associated Boundary object.
      *
      * \param boundary Boundary object.
      */
      void setBoundary(Boundary& boundary);

      /** 
      * Read parameters and initialize.
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);

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
  
      // Accessors

      /**
      * Return processor Grid by const reference.
      */
      const Grid& grid() const;

      /**
      * Get coordinate of processor within grid in specified direction.
      *
      * \param i index of Cartesian direction 0 <= i < Dimension
      */
      int gridCoordinate(int i) const;

      /**
      * Get dimension of processor grid in specified direction.
      *
      * \param i index of Cartesian direction 0 <= i < Dimension
      */
      int gridDimension(int i) const;

      /**
      * Is the grid periodic in the specified direction?
      *
      * \param i index of Cartesian direction 0 <= i < Dimension
      */
      bool gridIsPeriodic(int i) const;

      /**
      * Get rank of this processor in the processor grid.
      */
      int gridRank() const;

      /**
      * Is this the master processor (gridRank == 0) ?
      */
      bool isMaster() const;

      /**
      * Rank of the processor from which to receive for transfer (i, j).
      *
      * \param i index of Cartesian direction 0 <= i < Dimension
      * \param j index for positive (j = 0) or negative transfer (j = 1)
      */
      int sourceRank(int i, int j) const;

      /**
      * Rank of the processor to which to send for transfer (i, j).
      *
      * \param i index of Cartesian direction 0 <= i < Dimension
      * \param j index for positive (j = 0) or negative transfer (j = 1)
      */
      int destRank(int i, int j) const;

      /**
      * Shift to be applied to sent coordinates in transfer (i, j).
      *
      * When communicating an atomic position that is stored on the sending
      * processsor as x = {x[0], x[1], x[2]}, the value of coordinate x[i] 
      * should be sent as x[i] - shift(i, j)*length[i] during transfer (i, j), 
      * where length[i] is the length of the unit cell in direction i.
      *
      * \param i index of Cartesian direction 0 <= i < Dimension
      * \param j index for positive (j = 0) or negative (j = 1) transfer
      */
      int shift(int i, int j) const;

      #if UTIL_MPI
      /**
      * Return Cartesian communicator by reference.
      */
      MPI::Intracomm& communicator() const;
      #endif

      /**
      * Get one boundary of the domain owned by this processor.
      *
      * This processor owns a domain in which coordinate x[i] lies in 
      * the range domainBound(i, 0) < x[i] < domainBound(i, 1) for all
      * 0 <= i < Dimension.
      *
      * \param i index of Cartesian direction 0 <= i < Dimension
      * \param j index for lower (j=0) or upper (j=1) boundary.
      * \return boundary coordinate value
      */
      double domainBound(int i, int j) const;

      /**
      * Return rank of the processor whose domain contains a position.
      *
      * \param position position vector
      * \return rank of processor that owns this position.
      */
      int ownerRank(const Vector& position) const;

      /**
      * Is a position in the domain owned by this processor?
      *
      * \param position position vector
      * \return true if within the domain, false otherwise.
      */
      bool isInDomain(const Vector& position) const;

      /**
      * Has this Domain been initialized by calling readParam?
      */
      bool isInitialized() const;

      /**
      * Has an associated Boundary been set?
      */
      bool hasBoundary() const;

   private:

      // Ranks of sending (source) processor for each shift.
      FMatrix<int, Dimension, 2>  sourceRanks_;

      // Ranks of receiving (destination) processor for each shift.
      FMatrix<int, Dimension, 2>  destRanks_;

      // Shifts due to periodic boundary conditions in each send direction.
      FMatrix<int, Dimension, 2>  shift_;

      // Grid object with dimensions given by gridDimensions_.
      Grid grid_;

      // Number of processors in each direction.
      IntVector gridDimensions_;

      // Coordinates of this domain in processor grid.
      IntVector gridCoordinates_;

      // Rank of this processor in Cartesian communicator.
      int gridRank_;

      // Is each direction periodic (1 = true, 0 = false).
      FArray<bool, Dimension> gridIsPeriodic_;

      #if UTIL_MPI

      // Pointer to Intracommunicator.
      MPI::Intracomm* intracommPtr_;

      #endif

      // Pointer to associated Boundary object.
      Boundary* boundaryPtr_;

      // Is this object initialized (Has a grid been set?)
      bool isInitialized_;

      /*
      * Initialize after reading/loading gridDimensions.
      *
      * Used in implementation of readParameters and loadParameters.
      */
      void initialize();

   };

   #if UTIL_MPI

   /**
   * Return Cartesian communicator by reference.
   */
   inline MPI::Intracomm& Domain::communicator() const
   {
      assert(intracommPtr_);
      return *intracommPtr_; 
   }

   #endif

   /*
   * Return processor Grid by const reference.
   */
   inline const Grid& Domain::grid() const
   {  
      assert(isInitialized_);
      return grid_; 
   }

   /*
   * Return a coordinate of this processor within grid.
   */
   inline int Domain::gridCoordinate(int i) const
   {  
      assert(isInitialized_);
      return gridCoordinates_[i]; 
   }

   /**
   * Dimension of processor grid in specified direction.
   */
   inline int Domain::gridDimension(int i) const
   {  
      assert(isInitialized_);
      return gridDimensions_[i]; 
   }

   /**
   * Is the grid periodic in a specified direction?
   */
   inline bool Domain::gridIsPeriodic(int i) const
   {  
      assert(isInitialized_);
      return gridIsPeriodic_[i]; 
   }

   /*
   * Rank of this processor in Cartesian communicator.
   */
   inline int Domain::gridRank() const
   {  
      assert(isInitialized_);
      return gridRank_; 
   }

   /*
   * Is this the master processor (gridRank == 0)?
   */
   inline bool Domain::isMaster() const
   {  
      assert(isInitialized_);
      return bool(gridRank_ == 0); 
   }

   /*
   * Rank of the processor from which to receive in transfer (i, j).
   */
   inline int Domain::sourceRank(int i, int j) const
   {  
      assert(isInitialized_);
      return sourceRanks_(i, j); 
   }

   /*
   * Rank of the processor to which to send in transfer (i, j).
   */
   inline int Domain::destRank(int i, int j) const
   {  
      assert(isInitialized_);
      return destRanks_(i,j); 
   }

   /*
   * Shift applied in direction i for transfer (i, j)
   */
   inline int Domain::shift(int i, int j) const
   {  
      assert(isInitialized_);
      return shift_(i, j);  
   }

   /*
   * Has this Domain been initialized by calling readParam?
   */
   inline bool Domain::isInitialized() const
   {  return isInitialized_;  }

   /*
   * Has an associated Boundary been set ?
   */
   inline bool Domain::hasBoundary() const
   {  return (boundaryPtr_ != 0);  }

}
#endif
