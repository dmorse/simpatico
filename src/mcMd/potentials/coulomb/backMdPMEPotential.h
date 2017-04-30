#ifndef MD_PME_POTENTIAL_H
#define MD_PME_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/coulomb/MdCoulombPotential.h>            // base class
#include <mcMd/potentials/coulomb/EwaldRSpaceAccumulator.h>        // member
#include <mcMd/potentials/coulomb/EwaldInteraction.h>              // member

#include <util/space/IntVector.h>        // member template parameter
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
   * PME Coulomb potential class.
   *
   * This class implements the Particle Mesh for an Ewald  
   * summation of the Coulomb forces.
   *
   * \ingroup McMd_Coulomb_Module
   */
   class MdPMEPotential : public MdCoulombPotential
   {

   public:

      /**
      * Constructor.
      */
      MdPMEPotential(System& system);

      /**
      * Destructor (does nothing).
      */
      virtual ~MdPMEPotential();

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
      * Current number of wavevectors with |k| < kCutoff.
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
      * Compute kspace part of Coulomb pressure.
      *
      * \param stress (output) pressure
      */
      virtual void computeStress();
   
      //@}
      /// \name Accessors (const)
      //@{
      //
      EwaldRSpaceAccumulator& rSpaceAccumulator()
      {  return rSpaceAccumulator_; }

      EwaldInteraction& ewaldInteraction()
      {  return ewaldInteraction_; }

      //@}
      

   private:

      // 
      EwaldInteraction ewaldInteraction_;

      Simulation* simulationPtr_;

      System* systemPtr_;

      Boundary* boundaryPtr_;

      const Array<AtomType>* atomTypesPtr_;

      /// Grid Size.
      IntVector gridSize_;

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
      
      /// indicator.
      bool BCikinitialized_;

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

   private:

      // KSpace part of Coulomb energy
      Setable<double> kSpaceEnergy_;

      // KSpace part of Coulomb stress.
      Setable<Tensor> kSpaceStress_;

      /// Prefactor for self part coulomb energy.
      double selfPrefactor_;

      /// unitMatrix.
      Tensor unitTensor_;

   };


} 
#endif


