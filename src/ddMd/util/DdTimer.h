#ifndef DDMD_TIMER_H
#define DDMD_TIMER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>
#include <util/global.h>

namespace DdMd {

   /**
   * Class for measuring time intervals.
   *
   * Design adapted from the timer class in Lammps.
   */
   class DdTimer 
   {
   
   public:
   
      DdTimer(int size = 0);
      ~DdTimer();
      void clear();
      void start();
      void stamp(int id);
      void stop();
      double time(int id);
      double time();

      #ifdef UTIL_MPI
      void reduce(MPI::Intracomm& communicator);
      #endif
   
   private:
   
      Util::DArray<double> times_;
      double previous_;
      double begin_;
      double time_;
      int    size_;

   };

}
#endif
