#ifndef MCMD_MD_SPME_POTENTIAL_H
#define MCMD_MD_SPME_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/coulomb/MdCoulombPotential.h>      // base class
#include <mcMd/potentials/coulomb/EwaldRSpaceAccumulator.h>  // member
#include <simp/interaction/coulomb/EwaldInteraction.h>       // member

#include <util/space/IntVector.h>        // member template parameter
#include <util/space/Vector.h>           // member template parameter
#include <util/space/Tensor.h>           // member template parameter
#include <util/containers/Pair.h>        // member template parameter
#include <util/containers/GArray.h>      // member template
#include <util/containers/GridArray.h>   // member template
#include <util/misc/Setable.h>           // member template
#include <mcMd/chemistry/AtomType.h>     // member template parameter
#include <util/containers/Array.h>       // member class template
#include <util/boundary/Boundary.h>      // typedef

#include <complex>
#include <fftw3.h>

namespace McMd
{

   class Simulation;
   class System;

   typedef std::complex<double> DCMPLX;

   using namespace Util;
   using namespace Simp;

   /**
   * Smooth Particle-Mesh Ewald Coulomb potential for MD simulations.
   *
   * This class implements the smooth particle mesh ewald k-space
   * computations for the Coulomb energy and forces.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class MdSpmePotential : public MdCoulombPotential
   {

   public:

      /**
      * Constructor.
      */
      MdSpmePotential(System& system);

      /**
      * Destructor (destroy fftw plan).
      */
      virtual ~MdSpmePotential();

      /// \name Initialization
      //@{

      /**
      * Read parameters and initialize.
      *
      * \param in input stream
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

      //@}
      /// \name Interaction Parameters (get/set)
      //@{

      /**
      * Set a parameter value, identified by a string.
      *
      * \param name  parameter name
      * \param value  new value of parameter
      */
      void set(std::string name, double value);

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      */
      double get(std::string name) const;

      //@}
      /// \name System energy and stress.
      //@{

      /**
      * Precompute waves and influence function.
      */
      virtual void makeWaves();

      /**
      * Total number of waves or grid points.
      */
      virtual int nWave() const;

      /**
      * Add k-space Coulomb forces for all atoms.
      */
      virtual void addForces();

      /**
      * Calculate the long range kspace part of Coulomb energy.
      */
      virtual void computeEnergy();

      /**
      * place holder.
      */
      virtual void computeStress();

      //@}
      /// \name Miscellaneous Accessors
      //@{

      EwaldRSpaceAccumulator& rSpaceAccumulator()
      {  return rSpaceAccumulator_; }

      EwaldInteraction& ewaldInteraction()
      {  return ewaldInteraction_; }

      //@}

   private:

      // Ewald Interaction - core Ewald computations
      EwaldInteraction ewaldInteraction_;

      // Pointer to parent Simulation
      Simulation* simulationPtr_;

      // Pointer to parent System
      System* systemPtr_;

      // Pointer to boundary of associated System.
      Boundary* boundaryPtr_;

      // Pointer to array of atom types
      const Array<AtomType>* atomTypesPtr_;

      /// Grid dimensions - number of points in each direction
      IntVector gridDimensions_;

      /// Charge density assigned to r-space grid
      GridArray<DCMPLX> rhoR_;

      /// DFT of charge density on k-space grid
      GridArray<DCMPLX> rhoK_;

      /// Influence function
      GridArray<double> g_;
      
      /// Square magnitude of wavevectors
      GridArray<double> sqWaves_;
      
      /// Wavevectors
      GridArray<Vector> vecWaves_;
     
      /// Force grid x component
      GridArray<DCMPLX> xfield_;

      /// Force grid y component
      GridArray<DCMPLX> yfield_;

      /// Force grid z component 
      GridArray<DCMPLX> zfield_;

      /// order of basis spline
      int order_;
      
      /// FFT plan
      fftw_plan forward_plan;

      /// FFT plan for electric field
      fftw_plan xfield_backward_plan, yfield_backward_plan, zfield_backward_plan;

      /**
      * Set all elements of grid to 0.
      */
      template<class T> 
      void setGridToZero(GridArray<T>& grid);

      /**
      * Allocate and set all quantities that depend on grid dimensions.
      */
      void setGridDimensions();

      /**
      * Compute waves and influence function.
      */
      void computeWaves();

      /**
      * Compute components of b-factor in spme influence function.
      */
      double bfactor(double m , int dim);

      /**
      * Assign charges to grid points to compute rhoR_.
      */
      void assignCharges();

      /**
      * Expression for basis spline with order-5.
      */
      double basisSpline(double x);

   };

}
#endif
