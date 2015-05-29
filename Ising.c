//Compile: gcc Ising.c -lm -o Ising
//run: ./Ising


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int CalcE(int N,int arr[N][N]);
int Rel_CalcE(int N,int arr[N][N], int x, int y);
int flip(int N, int arr[N][N], float temp);
void main()
{
	FILE *fd;
	FILE *ft;
	fd=fopen("data.txt", "w");
	ft=fopen("timeC.txt", "w");
	clock_t start, end;
	double cpu_time_used;
	int N,t,i,j,num,E,k;
	float Eavg,temp;
	//N = 10;
	t = 100;

	for(N = 10; N<52;N+=2){
		printf("%d\n", N);
		int arr[N][N];
		Eavg = 0.;
		start = clock();
		for (i = 0; i<N; i++){
			for (j = 0; j<N; j++){
				num = rand() % 101;
				if(num<50){
					arr[i][j] = 1;
				}
				else{
					arr[i][j] = -1;
				}
			}

		}
		for(j = 0;j<51;j++){
			temp = j/10.;
			for(k = 0;k<N;k++){
				E = flip(N,arr, temp);
			}
			for (i = 0; i<t; i++){
				E = flip(N,arr, temp);
				Eavg +=E;
				//printf("%d\n", E);
			}
			Eavg = (Eavg/(N*N*t));
			printf("%f\n", Eavg);
			fprintf(fd, "%f\n", Eavg);
		}
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("%f\n", cpu_time_used);
		fprintf(ft, "%f\n", cpu_time_used);
		
	}
	fclose(fd);
	fclose(ft);


	//float msec = diff * 1000. / CLOCKS_PER_SEC;
	//printf("Time taken %f seconds", msec/1000.);
  	return;
}

int CalcE(int N,int arr[N][N])
{
	int M,i,j = 0;
	int E = 0;
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			M = M + arr[i][j];
			if(i<N-1){
				E+= -arr[i][j]*arr[i+1][j];
			}
			else if(i ==N-1){
				E+= -arr[i][j]*arr[0][j];
			}
			if(j<N-1){
				E+= -arr[i][j]*arr[i][j+1];
			}
			else if(j ==N-1){
				E+= -arr[i][j]*arr[i][0];
			}
		}
	}
	return E;
}

int Rel_CalcE(int N,int arr[N][N], int x, int y)
{
	//Calculate energy of nearest neighbors
	int E = 0;
	int i = x;
	int j = y;
	if(i>0 && i<N-1){
		E+= -arr[i][j]*arr[i+1][j]-arr[i][j]*arr[i-1][j];
	}
	else if(i ==N-1){
		E+= -arr[i][j]*arr[0][j]-arr[i][j]*arr[i-1][j];
	}
	else if(i ==0){
		E+= -arr[i][j]*arr[i+1][j]-arr[i][j]*arr[N-1][j];
	}
	if(j>0 && j<N-1){
		E+= -arr[i][j]*arr[i][j+1]-arr[i][j]*arr[i][j-1];
	}
	else if(j ==N-1){
		E+= -arr[i][j]*arr[i][0]-arr[i][j]*arr[i][j-1];
	}
	else if(j ==0){
		E+= -arr[i][j]*arr[i][j+1]-arr[i][j]*arr[i][N-1];
	}
	return E;
}

int flip(int N,int arr[N][N],float temp)
{
	int M,i,j,Enew,E0,E,x,y = 0;
	double p,num;
	int arrnew[N][N];
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			E0 = Rel_CalcE(N,arr,i,j);
			Enew = E0*-1;
			if(Enew<= E0){
				arr[i][j] *= -1;
			}
			else{
				p = exp(-(Enew-E0)/temp);
				num = (rand()% 101)/100.;
				if(p>=num){
					arr[i][j] *= -1;
				}
			}
		}
	}
	E = CalcE(N,arr);
	return E;
}


