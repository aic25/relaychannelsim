#include "const.h"
extern double CNR;
extern double channel_gain[N][2], p_in[KI][N];
double P_for_DD[KI][N][CL][CL][2], j[KI][CL][2];

extern errno_t err;	/* Error vector */

void RLS_initial_weigth(
	int receiver_index,
	double(*DFT_signal)[TS][N][2],
	double(*reference_signal)[KI][N][2],
	double(*W)[KI][CL][2],
	double(*H)[CL][BS][2]){

#ifndef RLS_VARIABLES
	int	n, np, v, ts, if_index = (receiver_index ? 0 : 1);
	double
		no,
		temp[2],
		temp2[2],
		temp3[2],
		power[2],
		gH[CL][2],
		power_i[2],
		R[CL][CL][2],
		e[N][NUMBER_OF_STEPS][2],
		tempvec[CL][2],
		tempvec2[CL][2],
		tempvec3[CL][2],
		tempmat[CL][CL][2],
		tempmat2[CL][CL][2];
#endif // !RLS_VARIABLES

#ifndef RLS_INITIALIZATION
	for (size_t f_index = 0; f_index < N; f_index++)
	{
		for (ts = 0; ts < CL; ts++){
			W[f_index][receiver_index][ts][0] = 0;
			W[f_index][receiver_index][ts][1] = 0;
			for (v = 0; v < CL; v++){
				if (ts == v)
				{
					P_for_DD[receiver_index][f_index][ts][v][0] = SIGMA_INVERSE;
				}
				else{
					P_for_DD[receiver_index][f_index][ts][v][0] = 0;
				}
				P_for_DD[receiver_index][f_index][ts][v][1] = 0;
				j[receiver_index][v][0] = 0;
				j[receiver_index][v][1] = 0;
			}
		}
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
				tempvec[ts][0] = DFT_signal[np][ts][f_index][0];
				tempvec[ts][1] = DFT_signal[np][ts][f_index][1];
			}			
			tempvec[CL - 1][0] = reference_signal[np][if_index][f_index][0];
			tempvec[CL - 1][1] = reference_signal[np][if_index][f_index][1];
			temp[0] = 0.0;
			temp[1] = 0.0;
			/*for (size_t cl = 0; cl < CL; cl++){
				temp[0] += W[f_index][0][cl][0] * DFT_signal[np][cl][f_index][0] - W[f_index][0][cl][1] * DFT_signal[np][cl][f_index][1];
				temp[1] += W[f_index][0][cl][0] * DFT_signal[np][cl][f_index][1] + W[f_index][0][cl][1] * DFT_signal[np][cl][f_index][0];
				}*/
			VectorMulVectorToScalar_CL(W[f_index][receiver_index], tempvec, temp);
			e[f_index][n][0] = reference_signal[np][receiver_index][f_index][0] - temp[0];
			e[f_index][n][1] = reference_signal[np][receiver_index][f_index][1] - temp[1];

			VectorConjugate_CL(tempvec, tempvec2);	/* y^C_m[n] */
			MatrixMulVectorToVector_CLxCL(P_for_DD[receiver_index][f_index], tempvec2, tempvec3);	/* P[n-1] x y^C_m[n] */
			VectorMulVectorToScalar_CL(tempvec, tempvec3, temp);	/* y^T_m[n] x P[n-1] x y^C_m[n] */
			temp[0] += FORGETTING_FACTOR;	/* lamda + y^T_m[n] x P[n-1] x y^C_m[n] */
			temp2[0] = 1.0; temp2[1] = 0.0;	/* {}^-1 */
			ScalarDivScalarToScalar(temp2, temp, temp3); /* {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */
			ScalarMulVectorToVector_CL(temp3, tempvec3, j[receiver_index]);	/*  y^H_m[n] x {lamda + y^T_m[n] x P[n-1] x y^C_m[n]}^-1 */

			temp2[0] = FF_INVERSE; temp2[1] = 0.0;	/* lamda^-1 */
			ScalarMulMatrixToMatrix_CLxCL(temp2, P_for_DD[receiver_index][f_index], tempmat);	/* lamda^-1 x P[n-1] */
			VectorMulMatrixToVector_CLxCL(tempvec, tempmat, tempvec2);	/* y^T_m[n] x lamda^-1 x P[n-1] */
			VectorMulVectorToMatrix_CLxCL(j[receiver_index], tempvec2, tempmat2);	/* g[n] x y^T_m[n] x lamda^-1 x P[n-1] */
			MatrixSubMatrixToMatrix_CL(tempmat, tempmat2, P_for_DD[receiver_index][f_index]);	/* lamda^-1 x P[n-1] - g[n] x y^T_m[n] x lamda^-1 x P[n-1] */

			ScalarMulVectorToVector_CL(e[f_index][n], j[receiver_index], tempvec);
			VectorAddVectorToVector_CL(W[f_index][receiver_index], tempvec, W[f_index][receiver_index]);

		}	/* end :: n */

#ifdef SD	 
		MatrixMulMatrixToMatrix_CLxCL(W[f_index], H[f_index], R);
		VectorConjugate_CL(W[f_index][receiver_index], gH);
		VectorMulVectorToScalar_CL(gH, W[f_index][receiver_index], power);
		power_i[0] = 0;
		for (size_t i = 1; i < M; i++){
			power_i[0] += sqr(R[0][i][0]) + sqr(R[0][i][1]);	/* 0.5=no */
		}
		channel_gain[f_index][0] = R[0][0][0];
		channel_gain[f_index][1] = -R[0][0][1];
		p_in[receiver_index][f_index] = (no)*power[0] + power_i[0];
#endif // SD


	}	/* end :: f_index */
}
