#include "GPUStructureFactorGrid.h"

#include "GPUStructureFactorGrid.cuh"
#include "cudacpu_vector_types.h"
#include "cudacpu_vector_functions.h"

namespace McMd
{

   using namespace Util;

   /// Constructor.
   GPUStructureFactorGrid::GPUStructureFactorGrid(System& system) 
    : StructureFactorGrid(system)
   {}


   /// Increment Structure Factor
   void GPUStructureFactorGrid::sample(long iStep)
   {
      if (isAtInterval(iStep))  {
  
         Vector  position;
         System::ConstMoleculeIterator  molIter;
         Molecule::ConstAtomIterator  atomIter;
         int  nSpecies, iSpecies, typeId, i, j;

         makeWaveVectors();

        
         // allocate temporary arrays
         int nType = system().simulation().nAtomType();

         float3 *h_wave_vectors = new float3[nWave_];
         float3 *h_pos = new float3[system().simulation().atomCapacity()];
         int *h_type = new int[system().simulation().atomCapacity()];
         float *h_mode = new float[nMode_*nType];
         float *h_sq = new float[nMode_*nWave_];

         // Load wave vectors
         for (int i = 0; i < nWave_; i++)
            h_wave_vectors[i] = make_float3((float) waveVectors_[i][0],
                                             (float) waveVectors_[i][1],
                                             (float) waveVectors_[i][2]);

         // Load atoms into temporary storage
         nSpecies = system().simulation().nSpecies();
         int nAtom = 0;
         for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
            system().begin(iSpecies, molIter);
            for ( ; molIter.notEnd(); ++molIter) {
               molIter->begin(atomIter);
               for ( ; atomIter.notEnd(); ++atomIter) {
                  position = atomIter->position();
                  typeId   = atomIter->typeId();

                  h_pos[nAtom] = make_float3(position[0],position[1],position[2]);
                  h_type[nAtom] = typeId;
                  nAtom++;
                  }

               }
            }

         // Load modes
         for (i = 0; i < nMode_; ++i) {
            for (j = 0; j < nType; ++j) {
               h_mode[nType*i + j] = modes_(i, j);
            }
         }
            
         // Calculate structure factors
         int res = gpu_sample_structure_factor(nWave_,
                                     h_wave_vectors,
                                     nAtom,
                                     h_pos,
                                     h_type,
                                     nType,
                                     nMode_,
                                     h_mode,
                                     h_sq,
                                     system().boundary().volume()
                                     );
          
         if (res)
             UTIL_THROW("Error updating structure factor.");

         // increment structure factors
         for (j = 0; j < nMode_; ++j) {
            double maxValue = 0.0;
            IntVector maxIntVector;
            double maxQ;
            for (i = 0; i < nWave_; ++i) {
               if ((double) h_sq[j*nWave_+i] >= maxValue) {
                  maxValue = (double) h_sq[j*nWave_+i];
                  maxIntVector = waveIntVectors_[i];
                  maxQ = waveVectors_[i].abs();
               }
               structureFactors_(i, j) += (double) h_sq[j*nWave_ + i];
            }
            maximumValue_[j].insert(maximumValue_[j].end(), 1, maxValue);
            maximumWaveIntVector_[j].insert(maximumWaveIntVector_[j].end(), 1, maxIntVector);
            maximumQ_[j].insert(maximumQ_[j].end(), 1, maxQ);
         }

         delete [] h_wave_vectors;
         delete [] h_pos;
         delete [] h_mode;
         delete [] h_type;
         delete [] h_sq;

          ++nSample_;

      }
   }    
}
