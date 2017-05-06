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
#include <mcMd/potentials/coulomb/EwaldInteraction.h>        // member

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

   /**
   * Ewald Coulomb potential class for MD simulations.
   *
   * This class implements the particle mesh ewald sum
   * of the Coulomb energy and forces.
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
      /// \name System energy and stress.
      //@{

      /**
      * place holder.
      */
      virtual void makeWaves();

      /**
      * place holder.
      */
      int nWave() const;

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
      { return ewaldInteraction_; }

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

      /// Grid Size.
      IntVector gridDimensions_;

      /// QGrid
      GridArray<DCMPLX> Qgrid_;

      /// Qhatgrid
      GridArray<DCMPLX> Qhatgrid_;

      /// BCGrid
      GridArray<double> BCgrid_;
      
      /// ik operator array. n-level rather than k-level ie. without prefactor 2Pi*I/L
      DArray<Vector> ikop_;

      /// force grid x component
      GridArray<DCMPLX> xfield_;

      /// force grid y component
      GridArray<DCMPLX> yfield_;

      /// force grid z component 
      GridArray<DCMPLX> zfield_;

      /// order of basis spline
      int order_;
      
      /// FFT plan.
      fftw_plan forward_plan;

      /// FFT plan for electric field.
      fftw_plan xfield_backward_plan, yfield_backward_plan, zfield_backward_plan;

      /**
       * set all elements to 0 for grid.
       */
      template<class T> 
      void initializeGrid(GridArray<T>& grid);

      /**
       * influence function, ie. BCgrid.
       */
      void influence_function();

      /**
       * compute components of B in BCgrid_.
       */
      double bfactor(double m , int dim);

      /**
       * ik operator.
       */
      void ik_differential_operator();
      
      /**
       * charge assignment function, ie. Qgrid_.
       */
      void spreadCharge();

      /**
       * expression for basis spline with order-5.
       */
      double basisSpline(double x);

      /// Prefactor for self-interaction correction.
      double selfPrefactor_;

   };

}
#endif
