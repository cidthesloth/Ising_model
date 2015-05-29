//-code=sm_20;
//Compile: nvcc Ising.cu -o Isingcuda
//Run: ./Isingcuda
#include <stdio.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <time.h>
#include <stdlib.h>

/*
#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)
*/

__global__ void init(unsigned int seed, curandState_t* states) {

  /* we have to initialize the state */
  curand_init(seed,threadIdx.x, 51,&states[threadIdx.x]);
}

__device__ int CalcE(int N,int arr[])
{        
       	int i,j;
	int M = 0;
	int E = 0;
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			M += arr[j+N*i];
			if(i<N-1){
				E+= -arr[j+N*i] * arr[j+1+N*i];
			}
			else if(i ==N-1){
				E+= -arr[j+N*i] * arr[j];
			}
			if(j<N-1){
				E+= -arr[j+N*i] * arr[j+1+N*i];
			}
			else if(j ==N-1){
				E+= -arr[j+N*i] * arr[j];
			}
		}
	}
	return E;
}

__device__ int Rel_CalcE(int N,int arr[] , int x, int y)
{
	//Calculate energy of nearest neighbors
	int E = 0;
	int i = x;
	int j = y;
	if(j>0 && j<N-1){
		E+= -arr[j+N*i] * arr[j+N*i+1]-arr[j+N*i] * arr[j+N*i-1];
	}
	else if(j ==N-1){
		E+= -arr[j+N*i] * arr[N*i]-arr[j+N*i] * arr[j+N*i-1];
	}
	else if(j ==0){
		E+= -arr[j+N*i] * arr[j+N*i+1]-arr[j+N*i] * arr[N*i+N-1];
	}
	if(i>0 && i<N-1){
		E+= -arr[j+N*i] * arr[j+N*(i+1)]-arr[j+N*i] * arr[j+N*(i-1)];
	}
	else if(i ==N-1){
		E+= -arr[j+N*i] * arr[i]-arr[j+N*i] * arr[j+N*(i-1)];
	}
	else if(i ==0){
		E+= -arr[j+N*i] * arr[j+N*(i+1)]-arr[j+N*i] * arr[j+N*(N-1)];
	}
	return E;

}

__device__ int flip(int N,int arr[],curandState_t* states)
{
	float p,num;
	float temp = threadIdx.x/10.;
	int i,j,Enew,E0,E= 0;
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			E0 = Rel_CalcE(N,arr,i,j);
			Enew = E0*(-1);
			if(Enew<= E0){
				arr[j+N*i] *= -1;
			}

			else{
				p = exp(-(Enew-E0)/temp);
				num = (curand(&states[threadIdx.x])% 101)/100.;
				if(p>=num){
					arr[j+N*i] *= -1;
				}			
			}				
		}
	}
	E = CalcE(N,arr);
	return E;
}

__global__ void Ising(int *arr, float *Earr,int *N,curandState_t* states)
{
  int itt = 100;
  int k = 0;
  int x = threadIdx.x;
  int y =(int ) *N;
  Earr[x]=0;

  for (k = 0; k<itt; k++){ 
  	Earr[x] += flip(y,arr);
  }
  Earr[x] = Earr[x]/itt;


}

int main()
{
    int i,j,num,k;
    int N;
    FILE *fd;
    FILE *ft;
    srand (time(NULL));
    curandState_t* states;
    fd=fopen("dataCuda.txt", "w");
    ft=fopen("timeCuda.txt", "w");
    for(N=10; N<52; N=N+2){

 	 /* allocate space on the GPU for the random states */
 	 cudaMalloc((void**) &states, N * sizeof(curandState_t));

  	 /* invoke the GPU to initialize all of the random states */
 	 init<<<1, 51>>>(time(0), states);

   	 size_t arrsize = N*N * sizeof(int);
   	 size_t Esize = 51 * sizeof(double);
	 size_t Nsize = sizeof(int);
   	 clock_t start, end;
   	 double cpu_time_used;
   	 start = clock();

   	 int* h_Arr = (int*)malloc(arrsize);
   	 float* h_E = (float*)malloc(51 * sizeof(float));


   	 // Allocate vectors in device memory
   	 int *d_Arr, *d_N;	
	 float *d_E;
   	 cudaMalloc(&d_Arr, arrsize);
   	 cudaMalloc(&d_E, Esize);
	 cudaMalloc(&d_N, Nsize);
	 //cudaCheckErrors("cudamalloc fail");

   	 for (i = 0; i<k; i++){
		for (j = 0; j<k; j++){
			num = rand() % 101;
			if(num<50){
				h_Arr[j+i*k] = 1;
			}
			else{
				h_Arr[j+i*k] = -1;
			}
			//printf("%d", h_Arr[j+i*k]);
			}

		}
   	 // Copy vectors from host memory to device memory
    	 cudaMemcpy(d_Arr, &h_Arr, arrsize, cudaMemcpyHostToDevice);
   	 cudaMemcpy(d_E, h_E, Esize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_N, &N, Nsize, cudaMemcpyHostToDevice);
	//cudaCheckErrors("cuda memcpy fail");

    	// Invoke kernel
   	 dim3 dimBlock( 51 );
   	 dim3 dimGrid( 1 );
    	Ising<<<dimGrid, dimBlock>>>(d_Arr, d_E, d_N, states);
   	 // Copy result from device memory to host memory
   	 cudaMemcpy(h_E, d_E, Esize, cudaMemcpyDeviceToHost);
    	cudaMemcpy(h_Arr, d_Arr, arrsize, cudaMemcpyDeviceToHost);
	//cudaCheckErrors("cudamemcpy or cuda kernel fail");
    	// Free device memory
   	cudaFree(d_Arr);
    	cudaFree(d_E);
	cudaFree(d_N);
	cudaFree(states);
    	end = clock();
    	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

	fprintf(fd, "%f\n", float(h_E[i]/(N*N*100.)));
	fprintf(ft, "%f\n", cpu_time_used);	
	}
	fclose(fd);
	fclose(ft);
    return EXIT_SUCCESS;
}
