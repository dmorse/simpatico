/* 
* Example of code added to a class header file to define a MPI type.
* Replace the string CPPTYPE with the appropriate class name.
*/


#ifdef UTIL_MPI
#include <util/mpi/MpiStructBuilder.h>
#endif

   #ifdef UTIL_MPI

      /**
      * Initialize MPI Datatype.
      */
      MPI::Datatype MpiTraits<CPPTYPE>::type = MPI::BYTE;
   
      /**
      * Commit MPI Datatype.
      */
      void CPPTYPE::commitMpiType() 
      {
         MpiStructBuilder builder;
         CPPTYPE           object;
   
         builder.setBase(&object);
         builder.addMember(&object.member1, MPI::DOUBLE);
          .....

         builder.commit(MpiTraits<CPPTYPE>::type);
      }

   #endif

