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
   * Plan for communication pattern for exchanging atoms
   * and communicating ghots.
   */
   class Plan
   {

   public:

      Plan() :
       flags_(0)
      {}

      void setExchange(int i, int j)
      { flags_ |= EMask[i][j]; }
 
      void clearExchange(int i, int j)
      { flags_ &= (~EMask[i][j]); }
 
      bool testExchange(int i, int j) const
      { return bool(flags_ & EMask[i][j]); }
 
      void setGhost(int i, int j)
      { flags_ |= GMask[i][j]; }
 
      void clearGhost(int i, int j)
      { flags_ &= (~GMask[i][j]); }
 
      bool testGhost(int i, int j) const
      { return bool(flags_ & GMask[i][j]); }
 
      unsigned int flags() const
      { return flags_; }
 
   private:

      unsigned int flags_;
 
      static unsigned int GMask[3][2];
      static unsigned int EMask[3][2];

   };


}
#endif
