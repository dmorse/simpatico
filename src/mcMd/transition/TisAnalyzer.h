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

      MdCriterion(MdSystem& system);

      setReactantCoordinate(double value);

      setProductCoordinate(double value);

      setInterfaceCoordinate(double value);

      /**
      * Return true iff simulation should stop.
      */
      virtual bool stopFlag() = 0;

      /**
      * Return true iff simulation should stop.
      */
      virtual bool saveFlag() = 0;

      /**
      * Return integer state code. 
      */
      virtual int stateId() = 0;

      /**
      * Return step index of last test.
      */
      int iStep();

   private:

      iStep_;

      stateId_;

      stopFlag_; 

      saveFlag_; 

      reactantCoordinate_;

      interfaceCoordinate_;

      productCoordinate_;

   };

}
#endif
