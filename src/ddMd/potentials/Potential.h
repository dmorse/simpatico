#ifndef DDMD_POTENTIAL_H
#define DDMD_POTENTIAL_H

#include <util/param/ParamComposite.h>  // base class
#include <util/misc/Setable.h>          // template for members
#include <util/space/Tensor.h>          // parameter for member

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   /**
   * A Potential represents a potential energy contribution.
   *
   * A Potential provides methods to calculate potential energy,
   * and associated atomic forces and stress.
   *  
   *  \ingroup DdMd_Potential_Module
   */
   class Potential : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Potential();

      /**
      * Destructor.
      */
      virtual ~Potential();

      /**
      * Set flag to identify if reverse communication is enabled.
      *
      * \param reverseUpdateFlag true if reverse communication is enabled.
      */
      void setReverseUpdateFlag(bool reverseUpdateFlag);

      /**
      * Get flag to identify if reverse communication is enabled.
      */
      bool reverseUpdateFlag() const;

      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add force contributions to all atomic forces.
      */
      virtual void computeForces() = 0;

      /**
      * Compute potential energy on all processors.
      *
      * This method must be called on all processors. The result is
      * stored on the master processor, and may be retrieved by 
      * calling energy() on this processor.
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator) = 0;
      #else
      virtual void computeEnergy() = 0;
      #endif

      /**
      * Return the total potential, from all processors.
      *
      * This method should only be called on the master (rank 0) 
      * processor, after a previous call to computeEnergy.
      */
      double energy() const;

      /**
      * Mark the energy as unknown (nullify).
      *
      * This should be called whenever atom positions or boundary changes.
      */
      void unsetEnergy();

      /**
      * Is the energy set (known)?
      */
      bool isEnergySet() const;

      /**
      * Compute stress on all processors.
      *
      * This method must be called on all processors. The result 
      * is stored on the master processor, and may be retrieved 
      * by calling stress() on this processor.
      */
      #ifdef UTIL_MPI
      virtual void computeStress(MPI::Intracomm& communicator)
      #else
      virtual void computeStress()
      #endif
      {}

      /**
      * Compute forces and stress for all processors.
      * 
      * Call on all processors. The default implementation just calls
      * computeForces() and computeStress() methods. Subclasses should
      * combine into a single loop.
      */
      #ifdef UTIL_MPI
      virtual void computeForcesAndStress(MPI::Intracomm& communicator);
      #else
      virtual void computeForcesAndStress();
      #endif

      /**
      * Return the stress tensor.
      *
      * This method should only be called on the master (rank 0) 
      * processor, after a previous call to computeStress.
      */
      Tensor stress() const;

      /**
      * Return the pressure. 
      *
      * This method should only be called on the master (rank 0) 
      * processor, after a previous call to computeStress.
      */
      double pressure() const;

      /**
      * Mark the stress as unknown (nullify).
      *
      * This should be called whenever atom positions or boundary changes.
      */
      void unsetStress();

      /**
      * Is the stress set (known)?
      */
      bool isStressSet() const;

      #ifdef UTIL_MPI
      /**
      * Is the potential in a valid internal state?
      *
      * Return true if valid, or throws an Exception.
      * Must be called on all processors.
      *
      * \param communicator domain communicator for all domain procs
      */
      virtual bool isValid(MPI::Intracomm& communicator) const;
      #endif

      //@}

   protected:

      /**
      * Set a value for the total energy.
      */
      void setEnergy(double energy);

      /**
      * Set a value for the total stress.
      */
      void setStress(const Tensor& stress);

      /**
      * Add a pair contribution to the stress tensor.
      *
      * \param f       pair force = f1 = -f2
      * \param dr      pair separation = r1 - r2 
      * \param stress  virial, incremented by dyad dr f
      */
      void incrementPairStress(const Vector& f, const Vector& dr, 
                               Tensor& stress) const;

      /**
      * Add local energies from all processors, set energy on master.
      *
      * Call on all processors.
      *
      * \param localEnergy  energy contribution from this processor
      * \param communicator domain communicator
      */
      #ifdef UTIL_MPI
      void reduceEnergy(double localEnergy, MPI::Intracomm& communicator);
      #else
      void reduceEnergy();
      #endif

      /**
      * Add local stresses from all processors, set total on master.
      *
      * Call on all processors.
      *
      * \param localStress  stress contribution from this processor
      * \param communicator domain communicator
      */
      #ifdef UTIL_MPI
      void reduceStress(Tensor& localStress, MPI::Intracomm& communicator);
      #else
      void reduceStress();
      #endif

   private:

      /// Total stress.
      Setable<Tensor> stress_;

      /// Total energy.
      Setable<double> energy_;

      /// Is reverse update communication enabled?
      bool reverseUpdateFlag_;

   };

   inline bool Potential::reverseUpdateFlag() const
   {  return reverseUpdateFlag_; }

   /*
   * Add a pair contribution to the virial tensor (protected).
   */
   inline void 
   Potential::incrementPairStress(const Vector& f, const Vector& dr, 
                                  Tensor& stress) const
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            stress(i, j) += f[i]*dr[j];
         }
      }
   }

}
#endif
