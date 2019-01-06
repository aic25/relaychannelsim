#include "const-function-declerations.h"		
extern double CNR, channel_power[N][2];
double channel_gain[N][2], p_in[N];

void MMSE_Weight_Matrix(double(*H)[TS][BS][2], double(*W)[BS][CL][2]){
	double no, hbhrH3[2], hbhbH3[2];
	int f_index;
	/* This one looks like numerical CNR calculation, while taking into account the coding rate */
	no = pow(10.0, -CNR / 10.0) / CODING_RATE;	/* Noise power. if you assume one for signal power... */

	for (f_index = 0; f_index < N; f_index++){
		ScalarMulScalar_H(H[f_index][0][2], H[f_index][1][2], hbhrH3);
		ScalarMulScalar_H(H[f_index][0][2], H[f_index][0][2], hbhbH3);
		hbhbH3[0] += no;
		ScalarDivScalarToScalar(hbhrH3, hbhbH3, W[f_index][0][0]);
	}
}