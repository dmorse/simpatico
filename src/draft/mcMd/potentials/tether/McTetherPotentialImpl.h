#ifdef  SIMP_TETHER
#ifndef MCMD_MC_TETHER_POTENTIAL_IMPL_H
#define MCMD_MC_TETHER_POTENTIAL_IMPL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/potentials/tether/McTetherPotential.h>       // base class
#include <mcMd/simulation/SubSystem.h>                      // base class
#include <mcMd/potentials/tether/McTetherPotentialImpl.h>
#include <util/global.h>

namespace Util
{
   class Vector;
   class Tensor;
}

namespace McMd
{

   using namespace Util;

   class System;

   /**
   * Template implementation of McTetherPotential.
   *
   * \ingroup McMd_Tether_Module
   */
   template <class Evaluator>
   class McTetherPotentialImpl : public McTetherPotential, public SubSystem
   {

   public:

      /** 
      * Constructor.
      */
      McTetherPotentialImpl(System& system);

      /** 
      * Constructor (copy from McTetherPotential)
      */
      McTetherPotentialImpl(McTetherPotentialImpl<Evaluator>& other);

      /** 
      * Destructor.
      */
      virtual ~McTetherPotentialImpl();

      /**
      * Read potential energy.
      * 
      * This method reads the tether potential Evaluator parameter
      * block. Before calling Evaluator::readParameters(), it passes
      * simulation().nTetherType() to Evaluator::setNAtomType().
      */
      virtual void readParameters(std::istream& in);

      /// \name Energy, Force, Stress Evaluators
      //@{

      /**
      * Returns potential energy for one tether.
      *
      * \param rSq  square of distance between tethered particles.
      * \param type type of tether.
      */
      virtual double energy(double rSq, int type) const;
  
      /**
      * Returns force/distance for one tether, for use in MD.
      *
      * A positive return value represents a repulsive radial force.
      *
      * \param rSq  square of distance between atom and anchor.
      * \param type type of tether.
      * \return     scalar repulsive force divided by distance.
      */
      virtual double forceOverR(double rSq, int type) const;
  
      /**
      * Return pair evaluator class name (e.g., "HarmonicTether").
      */
      virtual std::string evaluatorClassName() const;

      /**
      * Calculate the tether energy for one Atom.
      */
      void addForces();

      /**
      * Return total tether pair potential energy of this System.
      */
      double energy() const;

      //@}

      /**
      * Return tether evaluator by reference.
      */
      Evaluator& evaluator();

      /**
      * Return tether evaluator by const reference.
      */
      const Evaluator& evaluator() const;

   private:
  
      Evaluator* evaluatorPtr_;

      bool isCopy_;
 
      template <typename T>
      void computeStressImpl(T& stress) const;

   };

}

#include <mcMd/simulation/System.h> 
#include <mcMd/simulation/Simulation.h> 
#include <mcMd/species/Species.h>
#include <util/boundary/Boundary.h> 

#include <util/space/Dimension.h>
#include <util/space/Vector.h>

#include <fstream>

namespace McMd
{

   using namespace Util;

   /* 
   * Default constructor.
   */
   template <class Evaluator>
   McTetherPotentialImpl<Evaluator>::McTetherPotentialImpl(System& system)
    : McTetherPotential(),
      SubSystem(system),
      evaluatorPtr_(0),
      isCopy_(false)
   {  evaluatorPtr_ = new Evaluator(); }
 
   /* 
   * Constructor, copy from McTetherPotentialImpl<Evaluator>.
   */
   template <class Evaluator>
   McTetherPotentialImpl<Evaluator>::McTetherPotentialImpl(
                         McTetherPotentialImpl<Evaluator>& other)
    : McTetherPotential(),
      SubSystem(other.system())
      evaluatorPtr_(&other.evaluator()),
      isCopy_(true)
   {}
 
   /* 
   * Destructor. 
   */
   template <class Evaluator>
   McTetherPotentialImpl<Evaluator>::~McTetherPotentialImpl() 
   {
      if (evaluatorPtr_ && !isCopy_) {
         delete evaluatorPtr_;
      }
   }

   /* 
   * Read parameters from file.
   */
   template <class Evaluator>
   void McTetherPotentialImpl<Evaluator>::readParameters(std::istream &in) 
   {
      // Read only if not a copy.  Do not indent evaluator block.
      if (!isCopy_) {
         evaluator().setNAtomType(simulation().nAtomType());
         bool nextIndent = false;
         addParamComposite(evaluator(), nextIndent);
         evaluator().readParameters(in);
      }
   }
  
   /*
   * Return tether energy for a single tether.
   */
   template <class Evaluator>
   double McTetherPotentialImpl<Evaluator>::energy(double rSq, int typeId) 
      const
   { return evaluator().energy(rSq, typeId); }

   /*
   * Return force / separation for a single tether.
   */
   template <class Evaluator>
   double McTetherPotentialImpl<Evaluator>::forceOverR(double rSq, int typeId) const
   { return evaluator().forceOverR(double rSq, typeId); }

   /* 
   * Return total tether potential energy.
   */
   template <class Evaluator>
   double McTetherPotentialImpl<Evaluator>::energy() const
   { return 0.0; }

   /* 
   * Add tether forces to total.
   */
   template <class Evaluator>
   void McTetherPotentialImpl<Evaluator>::addForces()
   { }

   template <class Evaluator>
   inline Evaluator& McTetherPotentialImpl<Evaluator>::evaluator()
   { return *evaluatorPtr_; }

   template <class Evaluator>
   inline const Evaluator& McTetherPotentialImpl<Evaluator>::evaluator() const
   { return *evaluatorPtr_; }

   /*
   * Return tether potential evaluator class name.
   */
   template <class Evaluator>
   std::string McTetherPotentialImpl<Evaluator>::evaluatorClassName() const
   {  return evaluator().className(); }

}
#endif
#endif
