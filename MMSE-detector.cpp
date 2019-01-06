#include "const-function-declerations.h"	
extern double CNR, channel_gain[N][2], p_in[N], her[N][TS][2], channel_power[N][2], H[N][TS][BS][2];
extern int loop, R_index[N];

void CoherentDetection4SD(int data_number, double(*DFT_signal)[N][2], double(*WH)[BS][TS][2], double(*diLLR)[N_cbps]){
	int ts, n, index, i, offset;
	double temp[2], temp2[2], detected_signal[N_cbps], coef,
		no = pow(10.0, -CNR / 10.0) / CODING_RATE;

	index = 0;
	for (n = 0; n < N; n++) {
		temp[0] = temp[1] = 0.0;
		temp2[0] = temp2[1] = 0.0;

#ifndef MRC
		for (ts = 0; ts < TS; ts++) {

			temp[0] += WH[n][0][ts][0] * DFT_signal[ts][n][0] - WH[n][0][ts][1] * DFT_signal[ts][n][1];
			temp[1] += WH[n][0][ts][0] * DFT_signal[ts][n][1] + WH[n][0][ts][1] * DFT_signal[ts][n][0];
		}
#else 
		for (ts = 0; ts < TS; ts++) {	/* Recall that for MRC, optimal combining is the channel itself. :: her = h^H */
			temp[0] += her[n][ts][0] * DFT_signal[ts][n][0] - her[n][ts][1] * DFT_signal[ts][n][1];
			temp[1] += her[n][ts][0] * DFT_signal[ts][n][1] + her[n][ts][1] * DFT_signal[ts][n][0];
		}
#endif  		

		if (MODULATION == QPSK){
			/************************ QPSK ****************************/
			/*
			temp2[0] = (channel_gain[n][0] * temp[0] - channel_gain[n][1] * temp[1]) / p_in[n];
			temp2[1] = (channel_gain[n][0] * temp[1] + channel_gain[n][1] * temp[0]) / p_in[n];

			detected_signal[index] = temp2[0];
			detected_signal[index + 1] = temp2[1];
			*/
			detected_signal[index] = temp[0];
			detected_signal[index + 1] = temp[1];

			/***********************************************************/
		}
		else if (MODULATION == QAM16){
			coef = 2.0*OneBySqrt10 / (1.0 - channel_gain[n][0]);
			if (fabs(temp[0]) <= 2.0 * OneBySqrt10 * channel_gain[n][0]){
				detected_signal[index] = coef * OneBySqrt10 * temp[0];
			}
			else{
				if (temp[0] >= 0.0){
					detected_signal[index] = coef * OneBySqrt10 * (2.0 * temp[0] - 2.0 * OneBySqrt10 * channel_gain[n][0]);
				}
				else {
					detected_signal[index] = coef * OneBySqrt10 * (2.0 * temp[0] + 2.0 * OneBySqrt10 * channel_gain[n][0]);
				}
			}
			detected_signal[index + 1] = coef * OneBySqrt10 * (2.0 * OneBySqrt10 * channel_gain[n][0] - fabs(temp[0]));

			if (fabs(temp[1]) <= 2.0 * OneBySqrt10 * channel_gain[n][0]){
				detected_signal[index + 2] = coef * OneBySqrt10 * temp[1];
			}
			else{
				if (temp[1] >= 0.0){
					detected_signal[index + 2] = coef * OneBySqrt10 * (2.0 * temp[1] - 2.0 * OneBySqrt10 * channel_gain[n][0]);
				}
				else {
					detected_signal[index + 2] = coef * OneBySqrt10 * (2.0 * temp[1] + 2.0 * OneBySqrt10 * channel_gain[n][0]);
				}
			}
			detected_signal[index + 3] = coef * OneBySqrt10 * (2.0 * OneBySqrt10 * channel_gain[n][0] - fabs(temp[1]));


		}
		else if (MODULATION == QAM64){
			coef = 2.0*OneBySqrt42 / (1.0 - channel_gain[n][0]);
			for (i = 0; i < 2; i++){
				offset = i * 3;
				if (fabs(temp[i]) > 6.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset] = coef * (4.0 * temp[i] - 12.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						detected_signal[index + offset] = coef * (4.0 * temp[i] + 12.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else if (fabs(temp[i]) > 4.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset] = coef * (3.0 * temp[i] - 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						detected_signal[index + offset] = coef * (3.0 * temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else if (fabs(temp[i]) > 2.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset] = coef * (2.0 * temp[i] - 2.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						detected_signal[index + offset] = coef * (2.0 * temp[i] + 2.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else{
					detected_signal[index + offset] = coef * temp[i];
				}

				if (fabs(temp[i]) > 6.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 1] = coef * (-2.0 * temp[i] + 10.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						detected_signal[index + offset + 1] = coef * (2.0 * temp[i] + 10.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else if (fabs(temp[i]) > 2.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 1] = coef * (-temp[i] + 4.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else {
						detected_signal[index + offset + 1] = coef * (temp[i] + 4.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else{
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 1] = coef * (-2.0 * temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						detected_signal[index + offset + 1] = coef * (2.0 * temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}

				if (fabs(temp[i]) > 4.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 2] = coef * (-temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						detected_signal[index + offset + 2] = coef * (temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else{
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 2] = coef * (temp[i] - 2.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						detected_signal[index + offset + 2] = coef * (-temp[i] - 2.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
			}
		}
		index += MODULATION;

	}
	deinterleaver4SD(detected_signal, diLLR[data_number - Np]);
}

void CoherentDetection4HD(
	int data_number,
	double(*DFT_signal)[N][2],
	double(*W)[BS][CL][2],
	double(*H)[TS][BS][2],
	int(*diLLR)[Nd][N_cbps],
	int(*Decision),
	int(*R_index),
	int(*S_index)){

	int i, ts, n, index, detected_signal[KI][N_cbps], s3_index, s1_index, m1, m2, m3;
	double
		temp[2],
		temp1[2],
		e[2][2],
		vec_e[2][2],
		vec_eh[2][2],
		vec_u[2][2],
		y[TS][2],
		u[2],
		d[2],
		wu[2],
		S[3][2],
		herH[3][2],
		s_qpsk[4][2] = {
			{ OneBySqrt2, OneBySqrt2 },
			{ OneBySqrt2, -OneBySqrt2 },
			{ -OneBySqrt2, -OneBySqrt2 },
			{ -OneBySqrt2, OneBySqrt2 }
	};

	index = 0;

	for (n = 0; n < N; n++) {
		temp[0] = temp[1] = 0.0;
		temp1[0] = temp1[1] = 0.0;
		e[0][0] = 0.0; e[0][1] = 0.0;


		for (ts = 0; ts < TS; ts++){
			y[ts][0] = DFT_signal[ts][n][0];
			y[ts][1] = DFT_signal[ts][n][1];
		}
		VectorMulVectorToScalar_TS(her[n], y, d);
		VectorMulMatrixToVector_TSxCL(her[n], H[n], herH);

		if (MODULATION == QPSK){
			/************************ Hard QPSK decoding *************/
#ifdef MF_MLD
			for (m1 = 0; m1 < 4; m1++){
				for (m2 = 0; m2 < 4; m2++){
					for (m3 = 0; m3 < 4; m3++){
						S[0][0] = s_qpsk[m1][0];	S[0][1] = s_qpsk[m1][1];
						S[1][0] = s_qpsk[m2][0];	S[1][1] = s_qpsk[m2][1];
						S[2][0] = s_qpsk[m3][0];	S[2][1] = s_qpsk[m3][1];
						VectorMulVectorToScalar_CL(herH, S, u);
						e[1][0] = d[0] - u[0];
						e[1][1] = d[1] - u[1];
						if (m1 == 0 && m2 == 0 && m3 == 0){
							e[0][0] = e[1][0];
							e[0][1] = e[1][1];
							s1_index = 0;
							s3_index = 0;
						}
						else if (sqrt(sqr(e[1][0]) + sqr(e[1][1])) < sqrt(sqr(e[0][0]) + sqr(e[0][1]))){
							e[0][0] = e[1][0];
							e[0][1] = e[1][1];
							s1_index = m1;
							s3_index = m3;
						}
					}
				}
			}
#else
			for (m1 = 0; m1 < 4; m1++){
				for (m2 = 0; m2 < 4; m2++){
					for (m3 = 0; m3 < 4; m3++){
						S[0][0] = s_qpsk[m1][0];	S[0][1] = s_qpsk[m1][1];
						S[1][0] = s_qpsk[m2][0];	S[1][1] = s_qpsk[m2][1];
						S[2][0] = s_qpsk[m3][0];	S[2][1] = s_qpsk[m3][1];
						MatrixMulVectorToVector_3x2(H[n], S, vec_u);
						VectorSubVectorToVector_TS(y, vec_u, vec_e);
						//cout << "vec_e " << vec_e[0][0] << " " << vec_e[0][1] << endl;
						//cout << "vec_e " << vec_e[1][0] << " " << vec_e[1][1] << endl; getchar();
						//VectorConjugate_TS(vec_e, vec_eh);
						VectorMulVectorToScalarH_TS(vec_e, vec_e, e[1]);
						//cout << "e1 " << e[1][0] << " " << e[1][1] << endl; getchar();
						if (m1 == 0 && m2 == 0 && m3 == 0){
							e[0][0] = e[1][0];
							e[0][1] = e[1][1];
							s1_index = 0;
							s3_index = 0;
						}
						else if (sqrt(sqr(e[1][0]) + sqr(e[1][1])) < sqrt(sqr(e[0][0]) + sqr(e[0][1]))){
							e[0][0] = e[1][0];
							e[0][1] = e[1][1];
							s1_index = m1;
							s3_index = m3;
						}
					}
				}
			}
#endif // MF_MLD

			detected_signal[0][index] = (R_index[n] == 0) ? ((s_qpsk[s1_index][0] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][0] >= 0.0) ? 1 : 0);
			detected_signal[0][index + 1] = (R_index[n] == 0) ? ((s_qpsk[s1_index][1] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][1] >= 0.0) ? 1 : 0);

			detected_signal[1][index] = (R_index[n] == 1) ? ((s_qpsk[s1_index][0] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][0] >= 0.0) ? 1 : 0);
			detected_signal[1][index + 1] = (R_index[n] == 1) ? ((s_qpsk[s1_index][1] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][1] >= 0.0) ? 1 : 0);
			Decision[0] = (R_index[n] == 0) ? s1_index : s3_index; Decision[1] = (R_index[n] == 1) ? s1_index : s3_index;
			/************************ Hard QAM16 decoding *************/
		}
		else if (MODULATION == QAM16){
			detected_signal[R_index[n]][index] = (temp[0] >= 0.0) ? 1 : 0;
			detected_signal[R_index[n]][index + 1] = (-2.0*OneBySqrt10 <= temp[0] && temp[0] <= 2.0*OneBySqrt10) ? 1 : 0;
			detected_signal[R_index[n]][index + 2] = (temp[1] >= 0.0) ? 1 : 0;
			detected_signal[R_index[n]][index + 3] = (-2.0*OneBySqrt10 <= temp[1] && temp[1] <= 2.0*OneBySqrt10) ? 1 : 0;
		}

		/************************ Hard QAM64 Mod decoding *************/

		else if (MODULATION == QAM64){
			if (temp[0] > 0.0)detected_signal[R_index[n]][index] = 1;
			else detected_signal[R_index[n]][index] = 0;


			/*.....................*/
			if (fabs(temp[0]) > 4.0 / sqrt(42.0)){
				detected_signal[R_index[n]][index + 1] = 0;
				if (fabs(temp[0]) > 6.0 / sqrt(42.0)) detected_signal[R_index[n]][index + 2] = 0;
				else detected_signal[R_index[n]][index + 2] = 1;
			}
			/*.....................*/
			else{
				detected_signal[R_index[n]][index + 1] = 1;
				if (fabs(temp[0]) > 2.0 / sqrt(42.0)) detected_signal[R_index[n]][index + 2] = 1;
				else detected_signal[R_index[n]][index + 2] = 0;
			}

			/************************ ********************** *************/

			if (temp[1] > 0.0) detected_signal[R_index[n]][index + 3] = 1;
			else detected_signal[R_index[n]][index + 3] = 0;

			/*.....................*/
			if (fabs(temp[1]) > 4.0 / sqrt(42.0)){
				detected_signal[R_index[n]][index + 4] = 0;
				if (fabs(temp[1]) > 6.0 / sqrt(42.0)) detected_signal[R_index[n]][index + 5] = 0;
				else detected_signal[R_index[n]][index + 5] = 1;

			}
			/*.....................*/
			else{
				detected_signal[R_index[n]][index + 4] = 1;
				if (fabs(temp[1]) > 2.0 / sqrt(42.0)) detected_signal[R_index[n]][index + 5] = 1;
				else detected_signal[R_index[n]][index + 5] = 0;
			}
		}
		index += MODULATION;
	}
	if (CHANNEL_ENCODING) deinterleaver4HD(detected_signal[0], diLLR[0][data_number - Np]);
	else for (i = 0; i < N_dbps; i++) { diLLR[0][data_number - Np][i] = detected_signal[0][i]; }
	if (CHANNEL_ENCODING) deinterleaver4HD(detected_signal[1], diLLR[1][data_number - Np]);
	else for (i = 0; i < N_dbps; i++) { diLLR[1][data_number - Np][i] = detected_signal[1][i]; }
}
