#include "const-function-declerations.h"
extern double CNR;
extern double channel_gain[N][2], p_in[N], wh_sequence[BS][Np][N][2];
double P_for_DDWE[N][CL][CL][2], j[CL][2];

extern errno_t err;	/* Error vector */

void RLS_initial_weigth(
	int f_index,
	int receiver_index,
	double(*DFT_signal)[TS][N][2],
	double(*reference_signal)[TS][N][2],
	double(*W)[BS][CL][2],
	double(*H)[TS][BS][2],
	double(*e)[NUMBER_OF_STEPS][2]){

#ifndef RLS_VARIABLES
	int	n, np, v, ts;
	double
		sigma = 0.01,
		no,
		bir[2] = { 1, 0 },
		sifir[CL][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } },
		ff_inverse[2] = { FF_INVERSE, 0 },
		ff[2] = { FORGETTING_FACTOR, 0 },
		temp[2],
		temp2[2],
		temp3[2],
		power[2],
		power_i[2],
		gH[CL][2],
		x[CL][2],
		xh[CL][2],
		wh[CL][2],
		lpx[CL][2],
		tempvec[CL][2],
		tempvech[CL][2],
		tempvec2[CL][2],
		tempvec3[CL][2],
		R[BS][BS][2],
		lP_for_DDWE[N][CL][CL][2],
		tempmat[CL][CL][2],
		tempmat2[CL][CL][2];
#endif // !RLS_VARIABLES

#ifndef RLS_INITIALIZATION

	P_for_DDWE[f_index][0][0][0] = 1 / sigma; // / OneBySqrt2;
	P_for_DDWE[f_index][0][0][1] = 0;
	P_for_DDWE[f_index][0][1][0] = 0;
	P_for_DDWE[f_index][0][1][1] = 0;
	P_for_DDWE[f_index][0][2][0] = 0;
	P_for_DDWE[f_index][0][2][1] = 0;
	P_for_DDWE[f_index][1][0][0] = 0;
	P_for_DDWE[f_index][1][0][1] = 0;
	P_for_DDWE[f_index][1][1][0] = 1 / sigma;
	P_for_DDWE[f_index][1][1][1] = 0;
	P_for_DDWE[f_index][1][2][0] = 0;
	P_for_DDWE[f_index][1][2][1] = 0;
	P_for_DDWE[f_index][2][0][0] = 0;
	P_for_DDWE[f_index][2][0][1] = 0;
	P_for_DDWE[f_index][2][1][0] = 0;
	P_for_DDWE[f_index][2][1][1] = 0;
	P_for_DDWE[f_index][2][2][0] = 1 / sigma;//(DFT_signal[0][0][f_index][0] * DFT_signal[0][0][f_index][0] + DFT_signal[0][0][f_index][1] * DFT_signal[0][0][f_index][1]);
	P_for_DDWE[f_index][2][2][1] = 0;

	W[f_index][receiver_index][0][0] = 0;
	W[f_index][receiver_index][0][1] = 0;
	W[f_index][receiver_index][1][0] = 0;
	W[f_index][receiver_index][1][1] = 0;
	W[f_index][receiver_index][2][0] = 0;
	W[f_index][receiver_index][2][1] = 0;

	for (n = 0; n < Np + Nd; n++){	 // Np + Nd = Number of steps
		e[f_index][n][0] = 0;
		e[f_index][n][1] = 0;
	}
	temp[0] = 0.0;
	temp[1] = 0.0;
#endif // !RLS_INITIALIZATION

	/* -------------------------------------------- */
	no = pow(10.0, -CNR / 10.0) / CODING_RATE;	/* Noise power. Eb = 1 */

	for (n = 1; n < PREAMBLE_LENGTH + 1; n++){
		np = n - 1;
		for (ts = 0; ts < TS; ts++){
			x[ts][0] = DFT_signal[np][ts][f_index][0];
			x[ts][1] = DFT_signal[np][ts][f_index][1];
		}
		x[CL - 1][0] = wh_sequence[(receiver_index + 1) % 2][np][f_index][0];
		x[CL - 1][1] = wh_sequence[(receiver_index + 1) % 2][np][f_index][1];

		VectorConjugate_CL(x, xh);
		ScalarMulMatrixToMatrix_CLxCL(ff_inverse, P_for_DDWE[f_index], lP_for_DDWE[f_index]);
		MatrixMulVectorToVector_CLxCL(lP_for_DDWE[f_index], x, lpx);	// nom
		VectorMulVectorToScalar_CL(xh, lpx, temp2);
		temp2[0] += bir[0];	// denom
		ScalarDivScalarToScalar(bir, temp2, temp3);
		ScalarMulVectorToVector_CL(temp3, lpx, j);	// eq. 1

		VectorConjugate_CL(W[f_index][receiver_index], wh);
		VectorMulVectorToScalar_CL(x, wh, temp);
		
		e[f_index][np][0] = wh_sequence[receiver_index][np][f_index][0] - temp[0];
		e[f_index][np][1] = -(wh_sequence[receiver_index][np][f_index][1] - temp[1]);

		ScalarMulVectorToVector_CL(e[f_index][np], j, tempvec);
		VectorAddVectorToVector_CL(W[f_index][receiver_index], tempvec, tempvec2);
		VectorAddVectorToVector_CL(tempvec2, sifir, W[f_index][receiver_index]);

		VectorMulMatrixToVector_CLxCL(xh, lP_for_DDWE[f_index], tempvec);
		VectorMulVectorToMatrix_CLxCL(j, tempvec, tempmat);
		MatrixSubMatrixToMatrix_CL(lP_for_DDWE[f_index], tempmat, P_for_DDWE[f_index]);
	}	/* end :: n */


#ifdef SD
	MatrixMulMatrixToMatrix_BSxTS(W[f_index], H[f_index], R);
	VectorConjugate_TS(W[f_index][0], gH);
	VectorMulVectorToScalar_TS(gH, W[f_index][0], power);
	power_i[0] = 0;
	for (size_t i = 1; i < M; i++)
	{
		power_i[0] += sqr(R[0][i][0]) + sqr(R[0][i][1]);	/* 0.5=no */
	}
	channel_gain[f_index][0] = R[0][0][0];
	channel_gain[f_index][1] = -R[0][0][1];
	p_in[f_index] = (no)*power[0] + power_i[0];
#endif // SD


	/* end :: f_index */
}
