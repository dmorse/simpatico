#ifndef DDMD_TIMER_H
#define DDMD_TIMER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
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

      /**
      *  Clear all time statistics.
      */ 
      void clear();

      /**
      * Clear statistics, mark a start time.
      */ 
      void start();

      /**
      * Mark end of interval id.
      */ 
      void stamp(int id);

      /**
      * Stop total time accumulation.
      */ 
      void stop();

      /**
      * Get accumulated time for interval i, average per processor.
      */ 
      double time(int id) const;

      /**
      * Get total time since start time, average per processor.
      */ 
      double time() const;

      #ifdef UTIL_MPI
      /**
      * Upon return, times on every processor replaced by average over procs.
      */
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
