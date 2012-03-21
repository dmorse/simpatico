#ifndef MCMD_COMPOSITE_BOND_H
#define MCMD_COMPOSITE_BOND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/param/Begin.h>
#include <util/param/End.h>
#include <util/global.h>

#include <math.h>

namespace Util 
{  class Random; }

namespace McMd
{

   using namespace Util;

   /**
   * Composite bond potential.
   *
   * This class template creates a bond potential energy which is the 
   * sum of a BareBond bond potential and a BarePair pair potential.
   * This is intended for use in models in which the nonbonded pair 
   * interactions between bonded atoms are masked, to allow the pair
   * potential to be explicitly added back in.
   * 
   * \ingroup Potential_Module
   */
   template <class BareBond, class BarePair>
   class CompositeBond : public ParamComposite 
   {
   
   public:
   
      /**
      * Constructor.
      */
      CompositeBond();

      // Mutators

      /**
      * Read epsilon and sigma, initialize other variables.
      *
      * \pre nAtomType must be set, by calling setNAtomType().
      *
      * \param in  input stream 
      */
      void readParam(std::istream &in);
      
      /**  
      * Set nAtomType value.
      *
      * \param nBondType number of atom types.
      */
      void setNBondType(int nBondType);

      // Accessors
      
      /**
      * Returns interaction energy for a single pair of particles. 
      *
      * \param rsq    square of distance between pair of atoms
      * \param typeId type index for bond.
      * \return    pair interaction energy
      */
      double energy(double rsq, int typeId) const;
   
      /**
      * Returns ratio of bond force to atom separation.
      *
      * Multiply this quantity by the separation vector to obtain
      * the force vector. A positive value for the return value
      * represents a repulsive bond force.
      *
      * \param rsq        square of distance between particles
      * \param bondTypeId type index for bond.
      * \return force divided by distance 
      */
      double forceOverR(double rsq, int bondTypeId) const;

      /**
      * Throws exception if called.
      * 
      * \param randomPtr pointer to a random number generator
      * \param beta      inverse absolute temperature (inverse energy)
      * \param typeId    type index for bond.
      */  
      double randomBondLength(Random* randomPtr, double beta, int typeId) const;

      /**
      * Return name of instantiated class, with no spaces.
      */ 
      std::string className() const;
  
   private:

      // Bare bond potential object.
      BareBond bond_;

      // Bare pair potential object.
      BarePair pair_;
      
      // Class name string      
      std::string className_;
            
      /**
      * Copy constructor.
      */
      CompositeBond(const CompositeBond<BareBond, BarePair>& other);
 
   };
     
 
   /* 
   * Constructor.
   */
   template <class BareBond, class BarePair>
   CompositeBond<BareBond, BarePair>::CompositeBond() 
    : bond_(),
      pair_(),
      className_("CompositeBond")
   {
      className_ += "<";
      className_ += bond_.className();
      className_ += ",";
      className_ += pair_.className();
      className_ += ">";
   }
   
   
   /* 
   * Copy constructor.
   */   
   template <class BareBond, class BarePair>
   CompositeBond<BareBond, BarePair>::CompositeBond(const CompositeBond<BareBond, BarePair>& other) 
    : bond_(other.bond_),
      pair_(other.pair_)
   {}
     
   /* 
   * Read potential parameters from file.
   */  
   template <class BareBond, class BarePair>
   void CompositeBond<BareBond, BarePair>::readParam(std::istream &in)
   {
      std::string label;
      Begin* beginPtr = 0;
      End*   endPtr = 0;

      // Read bare bond potential
      beginPtr = &addBegin(bond_.className().c_str());
      beginPtr->setIndent(*this);
      beginPtr->readParam(in);
      readParamComposite(in, bond_);
      endPtr = &addEnd();
      endPtr->setIndent(*this);
      endPtr->readParam(in);

      // Read bare pair potential
      beginPtr = &addBegin(pair_.className().c_str());
      beginPtr->setIndent(*this);
      beginPtr->readParam(in);
      readParamComposite(in, pair_);
      endPtr = &addEnd();
      endPtr->setIndent(*this);
      endPtr->readParam(in);

   }
   
   /* 
   * Set nAtomType
   */
   template <class BareBond, class BarePair>
   void CompositeBond<BareBond, BarePair>::setNBondType(int nBondType) 
   {  
      bond_.setNBondType(nBondType);
      pair_.setNAtomType(1);
   }   

   /* 
   * Calculate interaction energy for a bond.
   */
   template <class BareBond, class BarePair>
   inline double 
   CompositeBond<BareBond, BarePair>::energy(double rsq, int typeId)
   const 
   {
      double total = bond_.energy(rsq, typeId);
      total += pair_.energy(rsq, 0, 0);
      return total;
   }
  
   /* 
   * Compute force/distance for a bond.
   */
   template <class BareBond, class BarePair>
   inline double 
   CompositeBond<BareBond, BarePair>::forceOverR(double rsq, int typeId)
   const 
   {
      double total = bond_.forceOverR(rsq, typeId);
      total += pair_.forceOverR(rsq, 0, 0);
      return total;
   }

   /*
   * Return name of instantiated class, with no spaces in template.
   */
   template <class BareBond, class BarePair>
   std::string CompositeBond<BareBond, BarePair>::className() const
   {  return className_; }
  
   /**
   * Throws exception if called.
   */  
   template <class BareBond, class BarePair>
   double CompositeBond<BareBond, BarePair>::randomBondLength(Random* randomPtr, double beta, int typeId) const
   {  UTIL_THROW("Unimplemented function");  }

}
#endif
