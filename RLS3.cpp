#include "const.h"
extern double CNR;
extern double channel_gain[N][2], p_in[N];
double P_for_DDWE[N][CL][CL][2], j[CL][2];

extern errno_t err;	/* Error vector */

void RLS_initial_weigth(
	int receiver_index,
	double(*DFT_signal)[TS][N][2],
	double(*reference_signal)[TS][N][2],
	double(*WH)[BS][CL][2],
	double(*H)[TS][BS][2],
	double(*wh_sequence)[Np][N],
	double(*e)[NUMBER_OF_STEPS][2]){

#ifndef RLS_VARIABLES
	int	n, np, v, ts;
	double
		sigma = 0.01,
		no,
		temp[2],
		gamma_inv[2],
		gamma[2],
		temp2[2],
		temp3[2],
		power[2],
		gH[CL][2],
		power_i[2],
		R[BS][BS][2],
		W[N][BS][CL][2],
		x[CL][2],
		xh[CL][2],
		k[CL][2],
		kh[CL][2],
		tempvec[CL][2],
		tempvec2[CL][2],
		tempvec3[CL][2],
		tempvec4[CL][2],
		tempmat[CL][CL][2],
		tempmat2[CL][CL][2],
		tempmat3[CL][CL][2];
#endif // !RLS_VARIABLES

#ifndef RLS_INITIALIZATION
	for (size_t f_index = 0; f_index < N; f_index++)
	{
		/*for (ts = 0; ts < CL; ts++){
			for (v = 0; v < CL; v++){
			if (ts == v){
			P_for_DDWE[f_index][ts][v][0] = SIGMA_INVERSE/OneBySqrt2;
			P_for_DDWE[f_index][ts][v][1] = SIGMA_INVERSE/OneBySqrt2;
			}
			else{
			P_for_DDWE[f_index][ts][v][0] = 0;
			P_for_DDWE[f_index][ts][v][1] = 0;
			}
			//P_for_DDWE[f_index][ts][v][1] = 0;
			j[v][0] = 0;
			j[v][1] = 0;
			}
			}  */

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
		WH[f_index][receiver_index][0][0] = 0;
		WH[f_index][receiver_index][0][1] = 0;
		WH[f_index][receiver_index][1][0] = 0;
		WH[f_index][receiver_index][1][1] = 0;
		WH[f_index][receiver_index][2][0] = 0;
		WH[f_index][receiver_index][2][1] = 0;
	}
	for (n = 0; n < Np + Nd; n++){
		for (size_t f_index = 0; f_index < N; f_index++){
			e[f_index][n][0] = 0;
			e[f_index][n][1] = 0;
		}
	}
	temp[0] = 0.0;
	temp[1] = 0.0;
#endif // !RLS_INITIALIZATION

	/* -------------------------------------------- */
	no = pow(10.0, -CNR / 10.0) / CODING_RATE;	/* Noise power. Eb = 1 */
	for (size_t f_index = 0; f_index < N; f_index++){
		for (n = 1; n < PREAMBLE_LENGTH + 1; n++){
			np = n - 1;
			for (ts = 0; ts < TS; ts++){
				x[ts][0] = DFT_signal[np][ts][f_index][0];
				x[ts][1] = DFT_signal[np][ts][f_index][1];
			}
			x[CL - 1][0] = wh_sequence[(receiver_index + 1) % 3][np][f_index];
			x[CL - 1][1] = 0;
			VectorConjugate_CL(x, xh);

			MatrixMulVectorToVector_CLxCL(P_for_DDWE[f_index], xh, tempvec);
			temp[0] = -FF_INVERSE; temp[1] = 0;
			ScalarMulVectorToVector_CL(temp, tempvec, k);	 // eq. 1		  
			VectorConjugate_CL(k, kh);


			VectorMulVectorToScalar_CL(k, x, temp);
			gamma_inv[0] = 1 - temp[0]; gamma_inv[1] = -temp[1]; // eq. 2
			temp[0] = 1; temp[1] = 0;
			ScalarDivScalarToScalar(temp, gamma_inv, gamma);

			VectorMulVectorToMatrix_CLxCL(kh, k, tempmat2);
			ScalarMulMatrixToMatrix_CLxCL(gamma, tempmat2, tempmat);
			temp[0] = FF_INVERSE; temp[1] = 0;
			ScalarMulMatrixToMatrix_CLxCL(temp, P_for_DDWE[f_index], tempmat2);
			MatrixSubMatrixToMatrix_CL(tempmat2, tempmat, P_for_DDWE[f_index]);	// eq. 3

			VectorMulVectorToScalar_CL(WH[f_index][receiver_index], x, temp);
			e[f_index][n][0] = wh_sequence[receiver_index][np][f_index] + temp[0];
			e[f_index][n][1] = 0 + temp[1];	// eq. 4

			temp[0] = e[f_index][n][0] * gamma[0] - e[f_index][n][1] * gamma[1];
			temp[1] = e[f_index][n][0] * gamma[1] + e[f_index][n][1] * gamma[0];	// eq. 5

			ScalarMulVectorToVector_CL(temp, k, tempvec);
			VectorAddVectorToVector_CL(WH[f_index][receiver_index], tempvec, tempvec2);	// eq. 6
			for (size_t cl = 0; cl < CL; cl++)
			{
				WH[f_index][receiver_index][cl][0] = tempvec2[cl][0];
				WH[f_index][receiver_index][cl][1] = tempvec2[cl][1];
			}
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


	}	/* end :: f_index */
}
