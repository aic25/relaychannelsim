#include "const-function-declerations.h"
extern double CNR;
extern double channel_gain[N][2], p_in[N], wh_sequence[BS][Np][N][2];
double P_for_ChDF[N][4][4][2], j_ChDF[4][2];
double sigma_old[N][2];

extern errno_t err;	/* Error vector */

void RLS_CombCoeff(
	int f_index,
	int receiver_index,
	int swap_index,
	double(*DFT_signal)[TS][N][2],
	double(*reference_signal)[TS][N][2],
	double(*W)[BS][CL][2],
	double(*H)[TS][BS][2],
	double(*e)[NUMBER_OF_STEPS][2],
	int First_Time){

#ifndef RLS_VARIABLES
	int	n, np, v, ts;
	double
		sigma = 0.001,
		no,
		bir[2] = { 1, 0 },
		sifir[CL][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } },
		ff_inverse[2] = { FF_INVERSE, 0 },
		ff[2] = { FORGETTING_FACTOR, 0 },
		temp[2],
		temp2[2],
		temp3[2],
		power[2],
		ffxs_old[2],
		power_i[2],
		gH[CL][2],
		y[CL][2],
		s_1[2],
		s_2[2],
		u[2],
		d[2],
		k[2],
		wu[2],
		ke[2],
		p_denom[2],
		p_nom[2],
		sigma_new[2];
	/*tempvec[CL][2],
	tempvech[CL][2],
	tempvec2[CL][2],
	tempvec3[CL][2],
	R[BS][BS][2],
	lP_for_DDWE[N][CL][CL][2],
	tempmat[CL][CL][2],
	tempmat2[CL][CL][2],   */
#endif // !RLS_VARIABLES

#ifndef RLS_INITIALIZATION
	if (First_Time) {
		sigma_old[f_index][0] = sigma;
		sigma_old[f_index][1] = 0;

		/*	P_for_DDWE[f_index][0][0][0] = 1 / sigma; // / OneBySqrt2;
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
			P_for_DDWE[f_index][2][2][1] = 0;		*/

		W[f_index][receiver_index][0][0] = 0;
		W[f_index][receiver_index][0][1] = 0;
		W[f_index][receiver_index][1][0] = 0;
		W[f_index][receiver_index][1][1] = 0;
		W[f_index][receiver_index][2][0] = 0;
		W[f_index][receiver_index][2][1] = 0;

		for (n = 0; n < Np + Nd; n++) {
			e[f_index][n][0] = 0;
			e[f_index][n][1] = 0;
		}
		temp[0] = 0.0;
		temp[1] = 0.0;
	}
#endif // !RLS_INITIALIZATION

	/* -------------------------------------------- */
	no = pow(10.0, -CNR / 10.0) / CODING_RATE;	/* Noise power. Eb = 1 */

	for (n = 1; n < PREAMBLE_LENGTH + 1; n++){
		np = n - 1;
		for (ts = 0; ts < TS; ts++){
			y[ts][0] = DFT_signal[np][ts][f_index][0];
			y[ts][1] = DFT_signal[np][ts][f_index][1];
		}
		if (swap_index == 0){				  
			u[0] = -(y[1][0] -
				(
				H[f_index][1][receiver_index][0] * wh_sequence[receiver_index][np][f_index][0] -
				H[f_index][1][receiver_index][1] * wh_sequence[receiver_index][np][f_index][1]
				) -
				(
				H[f_index][1][(receiver_index + 1) % 2][0] * wh_sequence[(receiver_index + 1) % 2][np][f_index][0] -
				H[f_index][1][(receiver_index + 1) % 2][1] * wh_sequence[(receiver_index + 1) % 2][np][f_index][1]
				));
			u[1] = -(y[1][1] -
				(
				H[f_index][1][receiver_index][0] * wh_sequence[receiver_index][np][f_index][1] +
				H[f_index][1][receiver_index][1] * wh_sequence[receiver_index][np][f_index][0]
				) -
				(
				H[f_index][1][(receiver_index + 1) % 2][0] * wh_sequence[(receiver_index + 1) % 2][np][f_index][1] +
				H[f_index][1][(receiver_index + 1) % 2][1] * wh_sequence[(receiver_index + 1) % 2][np][f_index][0]
				));
			if ((f_index == 5) && (DBG == ON))cout << "u: " << u[0] << " " << u[1] << endl;
			d[0] = y[0][0] -
				(
				H[f_index][0][receiver_index][0] * wh_sequence[receiver_index][np][f_index][0] -
				H[f_index][0][receiver_index][1] * wh_sequence[receiver_index][np][f_index][1]
				) -
				(
				H[f_index][0][(receiver_index + 1) % 2][0] * wh_sequence[(receiver_index + 1) % 2][np][f_index][0] -
				H[f_index][0][(receiver_index + 1) % 2][1] * wh_sequence[(receiver_index + 1) % 2][np][f_index][1]
				);
			d[1] = y[0][1] -
				(
				H[f_index][0][receiver_index][0] * wh_sequence[receiver_index][np][f_index][1] +
				H[f_index][0][receiver_index][1] * wh_sequence[receiver_index][np][f_index][0]
				) -
				(
				H[f_index][0][(receiver_index + 1) % 2][0] * wh_sequence[(receiver_index + 1) % 2][np][f_index][1] +
				H[f_index][0][(receiver_index + 1) % 2][1] * wh_sequence[(receiver_index + 1) % 2][np][f_index][0]
				);	
		}
		else{	 
			d[0] = y[1][0] -
				(
				H[f_index][1][receiver_index][0] * wh_sequence[receiver_index][np][f_index][0] -
				H[f_index][1][receiver_index][1] * wh_sequence[receiver_index][np][f_index][1]
				) -
				(
				H[f_index][1][(receiver_index + 1) % 2][0] * wh_sequence[(receiver_index + 1) % 2][np][f_index][0] -
				H[f_index][1][(receiver_index + 1) % 2][1] * wh_sequence[(receiver_index + 1) % 2][np][f_index][1]
				);
			d[1] = y[1][1] -
				(
				H[f_index][1][receiver_index][0] * wh_sequence[receiver_index][np][f_index][1] +
				H[f_index][1][receiver_index][1] * wh_sequence[receiver_index][np][f_index][0]
				) -
				(
				H[f_index][1][(receiver_index + 1) % 2][0] * wh_sequence[(receiver_index + 1) % 2][np][f_index][1] +
				H[f_index][1][(receiver_index + 1) % 2][1] * wh_sequence[(receiver_index + 1) % 2][np][f_index][0]
				);
			if ((f_index == 5) && (DBG == ON))cout << "u: " << u[0] << " " << u[1] << endl;
			u[0] = -(y[0][0] -
				(
				H[f_index][0][receiver_index][0] * wh_sequence[receiver_index][np][f_index][0] -
				H[f_index][0][receiver_index][1] * wh_sequence[receiver_index][np][f_index][1]
				) -
				(
				H[f_index][0][(receiver_index + 1) % 2][0] * wh_sequence[(receiver_index + 1) % 2][np][f_index][0] -
				H[f_index][0][(receiver_index + 1) % 2][1] * wh_sequence[(receiver_index + 1) % 2][np][f_index][1]
				));
			u[1] = -(y[0][1] -
				(
				H[f_index][0][receiver_index][0] * wh_sequence[receiver_index][np][f_index][1] +
				H[f_index][0][receiver_index][1] * wh_sequence[receiver_index][np][f_index][0]
				) -
				(
				H[f_index][0][(receiver_index + 1) % 2][0] * wh_sequence[(receiver_index + 1) % 2][np][f_index][1] +
				H[f_index][0][(receiver_index + 1) % 2][1] * wh_sequence[(receiver_index + 1) % 2][np][f_index][0]
				));
		}

		ScalarMulScalar_H(u, W[f_index][receiver_index][0], wu);
		e[f_index][np][0] = d[0] - wu[0];
		e[f_index][np][1] = d[1] - wu[1];

		if (First_Time) {

			ScalarMulScalar(ff, sigma_old[f_index], ffxs_old);
			sigma_new[0] = ffxs_old[0] + (sqr(u[0]) + sqr(u[1]));
			sigma_new[1] = ffxs_old[1];
			//if ((f_index == 5)&&(DBG==ON))cout << "div rls" << u[0] << "+j" << u[1] << " / " << sigma_new[0] << "+j" << sigma_new[1] << " = " << k[0] << "+j" << k[1] << endl; 
			ScalarDivScalarToScalar(u, sigma_new, k);

			ScalarMulScalar_H(k, e[f_index][np], ke);
			W[f_index][receiver_index][0][0] = W[f_index][receiver_index][0][0] + ke[0];
			W[f_index][receiver_index][0][1] = W[f_index][receiver_index][0][1] + ke[1];

			sigma_old[f_index][0] = sigma_new[0];
			sigma_old[f_index][1] = sigma_new[1];
		}
		//getchar();								  
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

