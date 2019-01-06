#include "const.h"
extern double CNR, EbNo;
extern double channel_gain[N][2], p_in[N], channel_power[N][2];
extern double mse_sum[LOOPN][NUMBER_OF_STEPS];

double P_for_DDWE2[N][FILTER_ORDER + 1][FILTER_ORDER + 1][2], j[FILTER_ORDER + 1][2];

extern int loop;
extern errno_t err;	/* Error vector */

void RLS_initial_weigth(int receiver_index, double(*DFT_signal)[TS][N][2], double(*reference_signal)[TS][N][2], double(*W)[BS][TS][2], double(*H)[TS][BS][2]){

#ifndef RLS_VARIABLES
	int	n, np, v, h, rx, f_index;
	double
		no,
		temp[2],
		temp2[2],
		temp3[2],
		power[2],
		gH[TS][2],
		power_i[2],
		R[BS][BS][2],
		e_stepav[N][2],
		mse[NUMBER_OF_STEPS],
		e[N][NUMBER_OF_STEPS][2],
		e_freqav[NUMBER_OF_STEPS][2],
		tempvec[FILTER_ORDER + 1][2],
		tempvec2[FILTER_ORDER + 1][2],
		tempvec3[FILTER_ORDER + 1][2],
		tempmat[FILTER_ORDER + 1][FILTER_ORDER + 1][2],
		tempmat2[FILTER_ORDER + 1][FILTER_ORDER + 1][2];
#endif // !RLS_VARIABLES

#ifndef RLS_INITIALIZATION
	for (size_t f_index = 0; f_index < N; f_index++)
	{
		for (h = 0; h < FILTER_ORDER + 1; h++){
			W[f_index][0][h][0] = 0;
			W[f_index][0][h][1] = 0;
			for (v = 0; v < FILTER_ORDER + 1; v++){
				if (h == v)
				{
					P_for_DDWE2[f_index][h][v][0] = SIGMA_INVERSE;
				}
				else{
					P_for_DDWE2[f_index][h][v][0] = 0;
				}
				P_for_DDWE2[f_index][h][v][1] = 0;
				j[v][0] = 0;
				j[v][1] = 0;
			}
		}
	}

	for (n = 0; n < Np + Nd; n++){
		e_freqav[n][0] = 0.0;
		e_freqav[n][1] = 0.0;
		for (f_index = 0; f_index < N; f_index++){
			e[f_index][n][0] = 0;
			e[f_index][n][1] = 0;
		}
	}

	for (f_index = 0; f_index < N; f_index++){
		e_stepav[f_index][0] = 0;
		e_stepav[f_index][1] = 0;
	}

	temp[0] = 0.0;
	temp[1] = 0.0;

	for (n = 0; n < Np + Nd; n++){
		mse[n] = 0.0;
		mse_sum[loop][n] = 0.0;
	}
#endif // !RLS_INITIALIZATION
	/* -------------------------------------------- */
	no = pow(10.0, -CNR / 10.0) / CODING_RATE;	/* Noise power. Eb = 1 */
	for (f_index = 0; f_index < N; f_index++){
		for (n = 1; n < PREAMBLE_LENGTH + 1; n++){
			np = n - 1;
			for (h = 0; h < TS; h++){
				tempvec[h][0] = DFT_signal[np][h][f_index][0];
				tempvec[h][1] = DFT_signal[np][h][f_index][1];
			}
			temp[0] = 0.0;
			temp[1] = 0.0;
			for (rx = 0; rx < TS; rx++){
				temp[0] += W[f_index][0][rx][0] * DFT_signal[np][rx][f_index][0] - W[f_index][0][rx][1] * DFT_signal[np][rx][f_index][1];
				temp[1] += W[f_index][0][rx][0] * DFT_signal[np][rx][f_index][1] + W[f_index][0][rx][1] * DFT_signal[np][rx][f_index][0];
			}
			e[f_index][n][0] = reference_signal[np][0][f_index][0] - temp[0];
			e[f_index][n][1] = reference_signal[np][0][f_index][1] - temp[1];

			VectorConjugate(tempvec, tempvec2);	/* y^C_m[n] */
			MatrixMulVectorToVector(P_for_DDWE2[f_index], tempvec2, tempvec3);	/* P[n-1] x y^C_m[n] */
			VectorMulVectorToScalar(tempvec, tempvec3, temp);	/* y^T_m[n] x P[n-1] x y^C_m[n] */
			temp[0] += FORGETTING_FACTOR;	/* lamda + y^T_m[n] x P[n-1] x y^C_m[n] */
			temp2[0] = 1.0; temp2[1] = 0.0;	/* {}^-1 */
			ScalarDivScalarToScalar(temp2, temp, temp3); /* {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */
			ScalarMulVectorToVector(temp3, tempvec3, j);	/*  y^H_m[n] x {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */

			temp2[0] = FF_INVERSE; temp2[1] = 0.0;	/* lamda^-1 */
			ScalarMulMatrixToMatrix(temp2, P_for_DDWE2[f_index], tempmat);	/* lamda^-1 x P[n-1] */
			VectorMulMatrixToVector(tempvec, tempmat, tempvec2);	/* y^T_m[n] x lamda^-1 x P[n-1] */
			VectorMulVectorToMatrix(j, tempvec2, tempmat2);	/* g[n] x y^T_m[n] x lamda^-1 x P[n-1] */
			MatrixSubMatrixToMatrix(tempmat, tempmat2, P_for_DDWE2[f_index]);	/* lamda^-1 x P[n-1] - g[n] x y^T_m[n] x lamda^-1 x P[n-1] */

			ScalarMulVectorToVector(e[f_index][n], j, tempvec);
			VectorAddVectorToVector(W[f_index][0], tempvec, W[f_index][0]);

		}	/* end :: n */

#ifdef SD
		MatrixMulMatrixToMatrix(W[f_index], H[f_index], R);
		VectorConjugate(W[f_index][0], gH);
		VectorMulVectorToScalar(gH, W[f_index][0], power);
		power_i[0] = sqr(R[0][1][0]) + sqr(R[0][1][1]);	/* 0.5=no */
		channel_gain[f_index][0] = R[0][0][0];
		channel_gain[f_index][1] = -R[0][0][1];
		p_in[f_index] = (no)*power[0] + power_i[0];
#endif // SD


	}	/* end :: f_index */
}
