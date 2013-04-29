#ifndef UTIL_MPI_LOADER_H
#define UTIL_MPI_LOADER_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/containers/DArray.h>      // template argument
#include <util/containers/FArray.h>      // template argument
#include <util/containers/DMatrix.h>     // template argument
#include <util/mpi/MpiFileIo.h>          // needed for template implementation
#include <util/mpi/MpiTraits.h>          // needed for template implementation
#include <util/mpi/MpiSendRecv.h>        // needed for template implementation
#include <util/global.h>


namespace Util
{

   /**
   * Methods to load data from input archive and broadcast to other processors.
   *
   * Each of the "load" method templates causes the ioProcessor to load a
   * variable from an archive and then broadcast it to all other processors
   * in the communicator provided by an associated MpiFileIo object.
   *
   *  \ingroup Mpi_Module
   */
   template <class IArchive>
   class MpiLoader 
   {

   public:

      /*
      * Constructor
      *  
      * \param mpiFileIo associated MpiFileIo object
      * \param archive   input archive from which data will be loaded
      */
      MpiLoader(MpiFileIo& mpiFileIo, IArchive& archive)
       : mpiFileIoPtr_(&mpiFileIo),
         archivePtr_(&archive) 
      {}

      /**  
      * Load and broadcast a single Type value.
      *
      * \param value  reference to a Type
      */
      template <typename Type>
      void load(Type &value);
   
      /**  
      * Load and broadcat a C array.
      *
      * \param value  pointer to array
      * \param n      number of elements
      */
      template <typename Type>
      void load(Type *value, int n);
   
      /**  
      * Load and broadcast a DArray < Type > container.
      *
      * \param array  DArray object
      * \param n      number of elements
      */
      template <typename Type>
      void load(DArray<Type>& array, int n);
   
      /**  
      * Load and broadcast an FArray < Type >.
      *
      * \param array  FArray object to be loaded
      */
      template <typename Type, int N>
      void load(FArray<Type, N >& array);
   
      /**  
      * Load and broadcast a 2D CArray of Type objects.
      *
      * \param value  pointer to array
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      */
      template <typename Type> void 
      load(Type *value, int m, int n);
  
      /**  
      * Load and broadcast a DMatrix<Type> object.
      *
      * \param matrix DMatrix object
      * \param m      number of rows (1st dimension)
      * \param n      number of columns (2nd dimension)
      */
      template <typename Type> 
      void load(DMatrix<Type>& matrix, int m, int n);

   private:

      // Pointer to associated MpiFileIo (passed to constructor).
      MpiFileIo*  mpiFileIoPtr_;

      // Pointer to associated input archive (passed to constructor).
      IArchive*  archivePtr_;
 
   };
 
   /*  
   * Load and broadcast a single Type value.
   */
   template <typename IArchive> 
   template <typename Type> 
   void MpiLoader<IArchive>::load(Type &value)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         *archivePtr_ >> value;
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Type>(mpiFileIoPtr_->ioCommunicator(), value, 0); 
      }
      #endif
   }
   
   /*  
   * Add a C array parameter, and load its elements.
   */
   template <typename IArchive> 
   template <typename Type>
   void MpiLoader<IArchive>::load(Type* value, int n)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         for (int i = 0; i < n; ++i) {
            *archivePtr_ >> value[i];
         }
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Type>(mpiFileIoPtr_->ioCommunicator(), value, n, 0); 
      }
      #endif
   }

   /*
   * Load a DArray < Type > container.
   */
   template <typename IArchive> 
   template <typename Type>
   void 
   MpiLoader<IArchive>::load(DArray<Type>& array, int n)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         *archivePtr_ >> array;
         if (array.capacity() < n) {
            UTIL_THROW("Error: DArray capacity < n");
         }
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Type>(mpiFileIoPtr_->ioCommunicator(), array, n, 0); 
      }
      #endif
   }

   /*
   * Load an FArray < Type, N > fixed-size array container.
   */
   template <typename IArchive> 
   template <typename Type, int N>
   void MpiLoader<IArchive>::load(FArray<Type, N >& array)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         for (int i = 0; i < N; ++i) {
            *archivePtr_ >> array[i];
         }
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Type>(mpiFileIoPtr_->ioCommunicator(), &(array[0]), N, 0); 
      }
      #endif
   }

   /*
   * Load a CArray2DParam < Type > C two-dimensional array parameter.
   */
   template <typename IArchive> 
   template <typename Type> 
   void 
   MpiLoader<IArchive>::load(Type *array, int m, int n)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         int i, j; 
         for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
               *archivePtr_ >> array[i*n + j];
            }
         }
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Type>(mpiFileIoPtr_->ioCommunicator(), &(array[0]), n*m, 0); 
      }
      #endif
   }
  
   /*
   * Add and load a DMatrix < Type > C two-dimensional matrix parameter.
   */
   template <typename IArchive> 
   template <typename Type> 
   void MpiLoader<IArchive>::load(DMatrix<Type>& matrix, int m, int n)
   {
      if (mpiFileIoPtr_->isIoProcessor()) {
         *archivePtr_ >> matrix;
      }
      #ifdef UTIL_MPI
      if (mpiFileIoPtr_->hasIoCommunicator()) {
         bcast<Type>(mpiFileIoPtr_->ioCommunicator(), matrix, m, n, 0); 
      }
      #endif
   }
  
} 
#endif
