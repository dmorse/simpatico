#include <cuda.h>
#include <stdio.h>
#include <cuComplex.h>

__global__ void kernel_calculate_sq_partial(
            int n_particles,
            cuComplex *fourier_mode_partial,
            float3 *pos,
            int n_wave,
            float3 *wave_vectors,
            float *d_modes, 
            int *d_type)
    {
    extern __shared__ cuComplex sdata[];

    unsigned int tidx = threadIdx.x;

    int j = blockIdx.x * blockDim.x + threadIdx.x;

    for (unsigned int i = 0; i < n_wave; i++) {
        float3 q = wave_vectors[i];

        cuComplex mySum = make_cuComplex(0.0,0.0);

        if (j < n_particles) {
            
            float3 p = pos[j];
            float dotproduct = q.x * p.x + q.y * p.y + q.z * p.z;
            int type = d_type[j];
            float mode = d_modes[type];
            cuComplex exponential = make_cuComplex(mode*cosf(dotproduct),
                                                   mode*sinf(dotproduct));
            mySum = cuCaddf(mySum,exponential);
        }
        sdata[tidx] = mySum;
    
       __syncthreads();

        // reduce in shared memory
        if (blockDim.x >= 512) {
           if (tidx < 256) {sdata[tidx] = mySum = cuCaddf(mySum,sdata[tidx+256]); }
            __syncthreads();
        }

        if (blockDim.x >= 256) {
           if (tidx < 128) {sdata[tidx] = mySum = cuCaddf(mySum, sdata[tidx + 128]); }
           __syncthreads(); 
        }

        if (blockDim.x >= 128) {
           if (tidx <  64) {sdata[tidx] = mySum = cuCaddf(mySum, sdata[tidx +  64]); }
           __syncthreads();
        }

        if (tidx < 32) {
            volatile cuComplex* smem = sdata;
            if (blockDim.x >= 64) {
                cuComplex rhs = cuCaddf(mySum, smem[tidx + 32]); 
                smem[tidx].x = rhs.x;
                smem[tidx].y = rhs.y;
                mySum = rhs;
            }
            if (blockDim.x >= 32) {
                cuComplex rhs = cuCaddf(mySum, smem[tidx + 16]); 
                smem[tidx].x = rhs.x;
                smem[tidx].y = rhs.y;
                mySum = rhs;
            }
            if (blockDim.x >= 16) {
                cuComplex rhs = cuCaddf(mySum, smem[tidx + 8]); 
                smem[tidx].x = rhs.x;
                smem[tidx].y = rhs.y;
                mySum = rhs;
            }
            if (blockDim.x >=  8) {
                cuComplex rhs = cuCaddf(mySum, smem[tidx + 4]); 
                smem[tidx].x = rhs.x;
                smem[tidx].y = rhs.y;
                mySum = rhs;
            }
            if (blockDim.x >=  4) {
                cuComplex rhs = cuCaddf(mySum, smem[tidx + 2]); 
                smem[tidx].x = rhs.x;
                smem[tidx].y = rhs.y;
                mySum = rhs;
            } 
            if (blockDim.x >=  2) {
                cuComplex rhs = cuCaddf(mySum, smem[tidx + 1]); 
                smem[tidx].x = rhs.x;
                smem[tidx].y = rhs.y;
                mySum = rhs;
            } 
        }

        // write result to global memeory
        if (tidx == 0)
           fourier_mode_partial[blockIdx.x + gridDim.x*i] = sdata[0];
    } // end loop over wave vectors
}

__global__ void kernel_calculate_norms(cuComplex* fourier_mode_partial,
                                       unsigned int nblocks, 
                                       float *sq_vec,
                                       int n_wave,
                                       float V)
    {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n_wave)
        return;

    // do final reduction of fourier mode
    cuComplex fourier_mode = make_cuComplex(0.0,0.0);
    for (unsigned int j = 0; j < nblocks; j++)
       fourier_mode = cuCaddf(fourier_mode, fourier_mode_partial[j + i*nblocks]);
 
    float normsq = fourier_mode.x * fourier_mode.x + fourier_mode.y * fourier_mode.y;
    sq_vec[i] = normsq/V;
    }

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
                                 ) 
    {
    float3* d_wave_vectors;
    float3* d_pos;
    int *d_type;
    float *d_modes;
    cuComplex *d_fourier_mode_partial;
    float *d_sq_vec;

    cudaError_t cudaStatus;

    cudaMalloc(&d_wave_vectors, sizeof(float3)*n_wave);
    cudaMemcpy(d_wave_vectors, h_wave_vectors, sizeof(float3)*n_wave, cudaMemcpyHostToDevice);

    cudaMalloc(&d_pos, sizeof(float3)*n_particles);
    cudaMemcpy(d_pos, h_pos, sizeof(float3)*n_particles, cudaMemcpyHostToDevice);

    cudaMalloc(&d_type, sizeof(int)*n_particles);
    cudaMemcpy(d_type, h_types, sizeof(int)*n_particles, cudaMemcpyHostToDevice);

    cudaMalloc(&d_modes, sizeof(float)*n_type*n_mode);
    cudaMemcpy(d_modes, h_modes, sizeof(float)*n_type*n_mode, cudaMemcpyHostToDevice);

    const unsigned int block_size_x = 256;
    unsigned int n_blocks_x = n_particles/block_size_x + 1;

    cudaMalloc(&d_fourier_mode_partial, sizeof(cuComplex)*n_wave*n_blocks_x);
    cudaMalloc(&d_sq_vec, sizeof(float)*n_wave*n_mode);

    for (int i = 0; i < n_mode; i++)
        {
        unsigned int shared_size = block_size_x * sizeof(cuComplex);
        kernel_calculate_sq_partial<<<n_blocks_x, block_size_x, shared_size>>>(
               n_particles,
               d_fourier_mode_partial,
               d_pos,
               n_wave,
               d_wave_vectors,
               d_modes + i*n_type,
               d_type);
 
        if (cudaStatus = cudaGetLastError()) {
               printf("CUDA ERROR (kernel_calculate_sq_partial): %s\n", cudaGetErrorString(cudaStatus));
               return 1;
        }

        // calculate final S(q) values of this mode
        const unsigned int block_size = 512;
        kernel_calculate_norms<<<n_wave/block_size + 1, block_size>>>(d_fourier_mode_partial,
                                                                      n_blocks_x,
                                                                      d_sq_vec + i*n_wave,
                                                                      n_wave,
                                                                      V);

        if (cudaStatus = cudaGetLastError())
            {
            printf("CUDA ERROR (kernel_calculate_norms): %s\n", cudaGetErrorString(cudaStatus));
            return 1;
            }


        } // end loop over modes

    // copy back structure factors
    cudaMemcpy(h_sq, d_sq_vec, n_wave*n_mode*sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_wave_vectors);
    cudaFree(d_pos);
    cudaFree(d_type);
    cudaFree(d_modes);
    cudaFree(d_fourier_mode_partial);
    cudaFree(d_sq_vec);

    return 0;
    }
