#ifndef MCMD_TIS_SNAPSHOT_ARRAY_H
#define MCMD_TIS_SNAPSHOT_ARRAY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

namespace McMd 
{

   class TisSnapShotArray
   {

      MdSnapShot(MdSystem&);

      reserve(int nStep);

      setNext(int iStep);

      clear();

      int nSnapShot();

      SnapShot& operator [] (int i);

   private:

      GArray<MdSnapShot> snapshots_;

   };

}
#endif
