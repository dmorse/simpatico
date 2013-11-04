#ifndef NVCC
#include "cudacpu_vector_types.h"
#endif

int gpu_sample_structure_factor(int n_wave,
                                 float3 *h_wave_vectors,
                                 unsigned int n_particles,
                                 float3 *h_pos,
                                 int *h_types,
                                 int n_type,   
                                 int n_mode,
                                 float *h_modes,
                                 float *h_sq,
                                 float V
                                 ); 

