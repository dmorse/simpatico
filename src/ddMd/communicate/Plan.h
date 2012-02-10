#ifndef PLAN_H
#define PLAN_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   /**
   * Plan for sending local or ghost atoms.
   */
   class Plan
   {

   public:

      Plan(){}

      void setExchange(int i, int j)
      { bits_ |= EMask[i][j]; }
 
      void clearExchange(int i, int j)
      { bits_ &= (~EMask[i][j]); }
 
      bool testExchange(int i, int j)
      { return bool(bits_ & ~EMask[i][j]); }
 
      void setGhost(int i, int j)
      { bits_ |= GMask[i][j]; }
 
      void clearGhost(int i, int j)
      { bits_ &= (~GMask[i][j]); }
 
      bool testGhost(int i, int j)
      { return bool(bits_ & ~GMask[i][j]); }
 
   private:

      unsigned int bits_;
 
      static unsigned int GMask[3][2];
      static unsigned int EMask[3][2];

   };


}
#endif
