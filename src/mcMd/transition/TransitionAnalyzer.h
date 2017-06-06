#ifndef MCMD_REACTION_ANALYZER_H
#define MCMD_REACTION_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd 
{

   class MdSystem;

   class ReactionAnalyzer : public SystemAnalyzer<MdSystem>
   {

      ReactionAnalyzer(MdSystem& system);

      void readParameter(std::istream& in);

      double computeReactionCoordinate()
      { return 0; }

      bool currentCoordinate() const
      {  return currentCoordinate_; }

      bool reactantCoordinate() const
      {  return reactantCoordinate_; }

      bool productCoordinate() const
      {  return productCoordinate_; }

      bool isReactant() const
      {  return currentCoordinate_ < reactantCoordinate_; }

      bool isProduct() const
      {  return currentCoordinate > productCoordinate_; }

      /**
      * Return step index of last test.
      */
      int iStep();
      {  return iStep_; }

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
      int iStep_;

   };

}
#endif
