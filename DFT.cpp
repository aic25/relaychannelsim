#include "const-function-declerations.h"

double IDFT_coef[N][N][2];


void DFT_initialize(){
	int n, k;
	double phase;

	for(k=0; k<N; k++) {		
		for(n=0; n<N; n++) {
			phase = (double) 2.0 * PI * n * k / N;
			IDFT_coef[n][k][0] = cos(phase);
			IDFT_coef[n][k][1] = sin(phase);
		}
	}
}


void IDFT(double (*input)[2], double (*output)[2]){
	int k, n;
	double temp[2];

	for(k=0; k<N; k++) {
		temp[0] = 0.0;
		temp[1] = 0.0;

		for(n=0; n<N; n++) {
			temp[0] += IDFT_coef[n][k][0] * input[n][0] - IDFT_coef[n][k][1] * input[n][1];
			temp[1] += IDFT_coef[n][k][0] * input[n][1] + IDFT_coef[n][k][1] * input[n][0];
		}
		output[k][0] = temp[0] / sqrt((double) N);
		output[k][1] = temp[1] / sqrt((double) N);
	}
}


void DFT(double (*input)[2], double (*output)[2]){
	int k, n;
	double temp[2];

	for(n=0; n<N; n++) {
		temp[0] = 0.0;
		temp[1] = 0.0;

		for(k=0; k<N; k++) {
			temp[0] += IDFT_coef[n][k][0] * input[k][0] + IDFT_coef[n][k][1] * input[k][1];
			temp[1] += IDFT_coef[n][k][0] * input[k][1] - IDFT_coef[n][k][1] * input[k][0];
		}
		output[n][0] = temp[0] / sqrt((double) N);
		output[n][1] = temp[1] / sqrt((double) N);
	}
}
