#ifndef MCMD_REACTION_ANALYZER_H
#define MCMD_REACTION_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/SystemAnalyzer.h>
#include <mcMd/mdSimulation/MdSystem.h>

namespace McMd 
{

   class MdSystem;

   class ReactionAnalyzer : public SystemAnalyzer<MdSystem>
   {

      /**
      * Constructor.
      *  
      * \param system associated MdSystem
      */
      ReactionAnalyzer(MdSystem& system);

      /**
      * Read parameters from file.
      * 
      * \param in  input file containing parameters
      */
      virtual void readParameter(std::istream& in);

      /**
      * Compute and return reaction coordinate value for associated MdSystem.
      */
      virtual double computeReactionCoordinate();

      /**
      * Sample state 
      * 
      * Compute and store reaction coordinate and store time stamp
      */
      virtual void sample(long iStep);

      /**
      * Get the most recently sampled value of reaction coordinate.
      */
      bool currentCoordinate() const;

      /**
      * Get maximum value of reaction coordinate for reactant domain.
      */
      bool reactantCoordinate() const;

      /**
      * Get minimum value of reaction coordinate for product domain.
      */
      bool productCoordinate() const;

      /**
      *
      * Is current reaction coordinate < reactant threshhold value?
      */
      bool isReactant() const;

      /**
      * Is current reaction coordinate > product threshhold value?
      */
      bool isProduct() const;

      /**
      * Return step index of most recent sample.
      */
      long iStep() const;

   protected:

      /** 
      * Value of reaction coordinate at most recent sample.
      */
      double currentCoordinate_;

      /**
      * Value of reaction coordinate below which system is reactant.
      */
      double reactantCoordinate_;

      /**
      * Value of reaction coordinate above which system is product.
      */
      double productCoordinate_;

      /**
      * Time step of last sample.
      */
      long iStep_;

   };

   inline
   bool ReactionAnalyzer::currentCoordinate() const
   {  return currentCoordinate_; }

   inline
   bool ReactionAnalyzer::reactantCoordinate() const
   {  return reactantCoordinate_; }

   inline
   bool ReactionAnalyzer::productCoordinate() const
   {  return productCoordinate_; }

   inline
   bool ReactionAnalyzer::isReactant() const
   {  return currentCoordinate_ < reactantCoordinate_; }

   inline
   bool ReactionAnalyzer::isProduct() const
   {  return currentCoordinate_ > productCoordinate_; }

   /*
   * Return step index of last test.
   */
   inline
   long ReactionAnalyzer::iStep() const
   {  return iStep_; }

}
#endif
