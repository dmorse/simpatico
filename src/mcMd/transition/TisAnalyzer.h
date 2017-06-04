#ifndef MCMD_TIS_ANALYZER_H
#define MCMD_TIS_ANALYZER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd 
{

   class MdSystem;

   class TisAnalyzer : public SystemAnalyzer<MdSystem>
   {

      TisAnalyzer(MdSystem& system);

      void readParameter(std::istream& in);

      void setLower(double lower);

      void setUpper(double lower);

      bool isReactant() const
      {  return isReactant_; }

      bool isProduct() const
      {  return isProduct_; }

      bool aboveLower() const
      {  return aboveLower_; }

      bool aboveUpper() const
      {  return aboveUpper_; }

      /**
      * Return step index of last test.
      */
      int iStep();
      {  return iStep_; }

   private:

      /** 
      * Value of reaction coordinate at most recent sample.
      */
      double currentRC_;

      /**
      * Value of reaction coordinate below which system is reactant.
      */
      double reactantRC_;

      /**
      * Value of reaction coordinate above which system is product.
      */
      double productRC_;

      /**
      * Lower threshold value for TIS algorithm.
      *
      * (All paths in ensemble must cross this).
      */
      double lowerRC_;

      /**
      * Upper threshold value for TIS algorithm.
      *
      * Paths terminate when they cross this.
      */
      double upperRC_;

      /**
      * Time step of last sample.
      */
      int iStep_;

      /**
      * Is the value in the reactant domain?
      */
      bool isReactant_;

      /**
      * Is the value in the product domain?
      */
      bool isProduct_;

      /**
      * Is the value above the lower TIS threshhold?
      */
      bool aboveLower_; 

      /**
      * Is the value above the upper TIS threshhold?
      */
      bool aboveUpper_; 

   };

}
#endif
