#ifndef DDMD_TIMER_CPP
#define DDMD_TIMER_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "DdTimer.h"

namespace DdMd
{

   DdTimer::DdTimer(int size)
   {
      times_.allocate(size);
      size_ = size;
      clear();
   }

   DdTimer::~DdTimer()
   {}

   void DdTimer::clear()
   {
      for (int i = 0; i < size_; i++) {
         times_[i] = 0.0;
      }
   }

   void DdTimer::start()
   {
      clear();
      begin_ = MPI_Wtime(); 
      previous_ = begin_;
   }

   void DdTimer::stamp(int id)
   {
      double current = MPI_Wtime();
      times_[id] += current - previous_;
      previous_ = current;
   }

   void DdTimer::stop()
   {  time_ = MPI_Wtime() - begin_; }

   #ifdef UTIL_MPI
   void DdTimer::reduce(MPI::Intracomm& communicator) 
   {
      int procs = communicator.Get_size();
      int rank  = communicator.Get_rank();
      bool isMaster = bool(rank == 0);
      double sum;
      for (int i = 0; i < size_; i++) {
         communicator.Reduce(&times_[i], &sum, 1, MPI::DOUBLE, MPI::SUM, 0);
         if (isMaster) times_[i] = sum/double(procs);
      }
      communicator.Reduce(&time_, &sum, 1, MPI::DOUBLE, MPI::SUM, 0);
      if (isMaster) time_ = sum/double(procs);
   }
   #endif

   double DdTimer::time(int id) const
   {  return times_[id]; }

   double DdTimer::time() const
   {  return time_; }

}
#endif
