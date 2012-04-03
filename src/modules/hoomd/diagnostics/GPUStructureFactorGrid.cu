#include <cuda.h>
#include <cublas_v2.h>
#include <stdio.h>

__global__ void kernel_load_mode_vec(
            int n_particles,
            float *d_modes, 
            int *d_type,
            cuComplex *d_mode_vec)
    {
    int n = blockIdx.x*blockDim.x + threadIdx.x;

    if (n >= n_particles)
        return;

    int type = d_type[n];
    d_mode_vec[n] = make_cuComplex(d_modes[type],0.0f);
    }

__global__ void kernel_load_matrix(
            int n_wave,
            cuComplex *exp_matrix,
            int pitch,
            float3 *pos,
            float3 *wave_vectors)
    {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= n_wave)
        return;

    int j = blockIdx.y;

    float3 q = wave_vectors[i];
    float3 p = pos[j];
    float dotproduct = q.x * p.x + q.y * p.y + q.z * p.z;
    // store in column-major format
    exp_matrix[j * pitch + i] = make_cuComplex(cosf(dotproduct),
                                                   sinf(dotproduct));
    }

__global__ void kernel_calculate_norms(cuComplex* fourier_mode_vec,
                                       float *sq_vec,
                                       int n_wave,
                                       float V)
    {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= n_wave)
        return;
    cuComplex fourier_mode = fourier_mode_vec[i];
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
    cuComplex *d_mode_vec;
    cuComplex *d_fourier_mode_vec;
    float *d_sq_vec;

    cudaError_t cudaStatus;

    cudaMalloc(&d_wave_vectors, sizeof(float3)*n_wave);
    cudaMemcpy(d_wave_vectors, h_wave_vectors, sizeof(float3)*n_wave, cudaMemcpyHostToDevice);

    cudaMalloc(&d_pos, sizeof(float3)*n_particles);
    cudaMemcpy(d_pos, h_pos, sizeof(float3)*n_particles, cudaMemcpyHostToDevice);

    cudaMalloc(&d_type, sizeof(int)*n_particles);
    cudaMemcpy(d_type, h_types, sizeof(int)*n_particles, cudaMemcpyHostToDevice);

    cudaMalloc(&d_modes, sizeof(float)*n_type*n_mode);
    cudaMemcpy(d_modes, h_modes,sizeof(float)*n_type*n_mode,cudaMemcpyHostToDevice);

    cudaMalloc(&d_mode_vec, sizeof(cuComplex)*n_particles);
    cuComplex *d_exp_matrix;
    size_t pitch;

    cudaMallocPitch((void **)&d_exp_matrix, &pitch, (size_t) (sizeof(cuComplex)*n_wave), (size_t)n_particles);
    pitch/=sizeof(cuComplex);
    cudaMalloc(&d_fourier_mode_vec, sizeof(cuComplex)*n_wave);
    cudaMalloc(&d_sq_vec, sizeof(float)*n_wave*n_mode);


    // initialize cuBLAS
    cublasStatus_t stat;
    cublasHandle_t handle;
    stat = cublasCreate(&handle);
    if ( stat != CUBLAS_STATUS_SUCCESS )
        {
        printf("CUBLAS Error %d\n", stat);
        return 1;
        }

    for (int i = 0; i < n_mode; i++)
        {
        // load mode vector
        int block_size = 512;

        kernel_load_mode_vec<<<n_particles/block_size + 1, block_size>>>(
            n_particles,
            d_modes + i*n_type, 
            d_type,
            d_mode_vec);

        if (cudaStatus = cudaGetLastError())
            {
            printf("CUDA ERROR: %s\n", cudaGetErrorString(cudaStatus));
            return 1;
            }

        // load exponential factor matrix
        dim3 dimGrid(n_wave/block_size + 1,n_particles);
        dim3 dimBlock(block_size,1);
        kernel_load_matrix<<<dimGrid, dimBlock>>>(
            n_wave,
            d_exp_matrix,
            pitch,
            d_pos,
            d_wave_vectors);

        if (cudaStatus = cudaGetLastError())
            {
            printf("CUDA ERROR: %s\n", cudaGetErrorString(cudaStatus));
            return 1;
            }

        // matrix multiplication of exp_matrix with mode_vec
        cuComplex alpha = make_cuComplex(1.0f,0.0f);

        cuComplex beta = make_cuComplex(0.0f,0.0f);

        stat = cublasCgemv(handle,
                    CUBLAS_OP_N,
                    n_wave, 
                    n_particles,
                    &alpha,
                    d_exp_matrix,
                    pitch,
                    d_mode_vec, 1,
                    &beta,
                    d_fourier_mode_vec, 1);

        if ( stat != CUBLAS_STATUS_SUCCESS )
            {
            printf("CUBLAS Error %d\n", stat);
            return 1;
            }

        // calculate norms of the entries of sq_vec
        
        kernel_calculate_norms<<<n_wave/block_size + 1, block_size>>>(d_fourier_mode_vec,
                                                                      d_sq_vec + i*n_wave,
                                                                      n_wave,
                                                                      V);

        if (cudaStatus = cudaGetLastError())
            {
            printf("CUDA ERROR: %s\n", cudaGetErrorString(cudaStatus));
            return 1;
            }


        } // end loop over modes

    stat = cublasDestroy ( handle ) ;
    if ( stat != CUBLAS_STATUS_SUCCESS )
        {
        printf("CUBLAS Error %d\n", stat);
        return 1;
        }
    // copy back structure factors
    cudaMemcpy(h_sq, d_sq_vec, n_wave*n_mode*sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(d_wave_vectors);
    cudaFree(d_pos);
    cudaFree(d_type);
    cudaFree(d_modes);
    cudaFree(d_mode_vec);
    cudaFree(d_exp_matrix);
    cudaFree(d_fourier_mode_vec);
    cudaFree(d_sq_vec);

    return 0;
    }
