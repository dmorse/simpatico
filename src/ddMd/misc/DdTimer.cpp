/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
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
      time_ = 0.0;
   }

   void DdTimer::start()
   {
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
   {  time_ += MPI_Wtime() - begin_; }

   #ifdef UTIL_MPI
   void DdTimer::reduce(MPI_Comm communicator) 
   {
      int procs;
      MPI_Comm_size(communicator, &procs);
      double sum;
      for (int i = 0; i < size_; i++) {
         MPI_Allreduce(&times_[i], &sum, 1, MPI_DOUBLE, MPI_SUM, communicator);
         times_[i] = sum/double(procs);
      }
      MPI_Allreduce(&time_, &sum, 1, MPI_DOUBLE, MPI_SUM, communicator);
      time_ = sum/double(procs);
   }
   #endif

   double DdTimer::time(int id) const
   {  return times_[id]; }

   double DdTimer::time() const
   {  return time_; }

}
