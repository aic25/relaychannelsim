#include "const.h"
extern double channel_gain[N][2], p_in[N], her[N][TS][2], channel_power[N][2];
extern int loop;

void CoherentDetection4SD(int data_number, double(*DFT_signal)[N][2], double(*WH)[BS][TS][2], double(*diLLR)[N_cbps]){
	int ts, n, index, i, offset;
	double temp[2], temp2[2], Detected_Signal[N_cbps], coef;

	index = 0;
	for (n = 0; n < N; n++) {
		temp[0] = temp[1] = 0.0;
		temp2[0] = temp2[1] = 0.0;

		for (ts = 0; ts < TS; ts++) {
			temp[0] += WH[n][0][ts][0] * DFT_signal[ts][n][0] - WH[n][0][ts][1] * DFT_signal[ts][n][1];
			temp[1] += WH[n][0][ts][0] * DFT_signal[ts][n][1] + WH[n][0][ts][1] * DFT_signal[ts][n][0];
		}


		if (MODULATION == QPSK){
			/************************ QPSK ****************************/
			temp2[0] = (channel_gain[n][0] * temp[0] - channel_gain[n][1] * temp[1]) / p_in[n];
			temp2[1] = (channel_gain[n][0] * temp[1] + channel_gain[n][1] * temp[0]) / p_in[n];

			Detected_Signal[index] = temp2[0];
			Detected_Signal[index + 1] = temp2[1];
			/***********************************************************/
		}
		else if (MODULATION == QAM16){
			coef = 2.0*OneBySqrt10 / (1.0 - channel_gain[n][0]);
			if (fabs(temp[0]) <= 2.0 * OneBySqrt10 * channel_gain[n][0]){
				Detected_Signal[index] = coef * OneBySqrt10 * temp[0];
			}
			else{
				if (temp[0] >= 0.0){
					Detected_Signal[index] = coef * OneBySqrt10 * (2.0 * temp[0] - 2.0 * OneBySqrt10 * channel_gain[n][0]);
				}
				else {
					Detected_Signal[index] = coef * OneBySqrt10 * (2.0 * temp[0] + 2.0 * OneBySqrt10 * channel_gain[n][0]);
				}
			}
			Detected_Signal[index + 1] = coef * OneBySqrt10 * (2.0 * OneBySqrt10 * channel_gain[n][0] - fabs(temp[0]));

			if (fabs(temp[1]) <= 2.0 * OneBySqrt10 * channel_gain[n][0]){
				Detected_Signal[index + 2] = coef * OneBySqrt10 * temp[1];
			}
			else{
				if (temp[1] >= 0.0){
					Detected_Signal[index + 2] = coef * OneBySqrt10 * (2.0 * temp[1] - 2.0 * OneBySqrt10 * channel_gain[n][0]);
				}
				else {
					Detected_Signal[index + 2] = coef * OneBySqrt10 * (2.0 * temp[1] + 2.0 * OneBySqrt10 * channel_gain[n][0]);
				}
			}
			Detected_Signal[index + 3] = coef * OneBySqrt10 * (2.0 * OneBySqrt10 * channel_gain[n][0] - fabs(temp[1]));


		}
		else if (MODULATION == QAM64){
			coef = 2.0*OneBySqrt42 / (1.0 - channel_gain[n][0]);
			for (i = 0; i < 2; i++){
				offset = i * 3;
				if (fabs(temp[i]) > 6.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						Detected_Signal[index + offset] = coef * (4.0 * temp[i] - 12.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						Detected_Signal[index + offset] = coef * (4.0 * temp[i] + 12.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else if (fabs(temp[i]) > 4.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						Detected_Signal[index + offset] = coef * (3.0 * temp[i] - 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						Detected_Signal[index + offset] = coef * (3.0 * temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else if (fabs(temp[i]) > 2.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						Detected_Signal[index + offset] = coef * (2.0 * temp[i] - 2.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						Detected_Signal[index + offset] = coef * (2.0 * temp[i] + 2.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else{
					Detected_Signal[index + offset] = coef * temp[i];
				}

				if (fabs(temp[i]) > 6.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						Detected_Signal[index + offset + 1] = coef * (-2.0 * temp[i] + 10.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						Detected_Signal[index + offset + 1] = coef * (2.0 * temp[i] + 10.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else if (fabs(temp[i]) > 2.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						Detected_Signal[index + offset + 1] = coef * (-temp[i] + 4.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else {
						Detected_Signal[index + offset + 1] = coef * (temp[i] + 4.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else{
					if (temp[i] >= 0.0){
						Detected_Signal[index + offset + 1] = coef * (-2.0 * temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						Detected_Signal[index + offset + 1] = coef * (2.0 * temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}

				if (fabs(temp[i]) > 4.0 * OneBySqrt42 * channel_gain[n][0]){
					if (temp[i] >= 0.0){
						Detected_Signal[index + offset + 2] = coef * (-temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						Detected_Signal[index + offset + 2] = coef * (temp[i] + 6.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
				else{
					if (temp[i] >= 0.0){
						Detected_Signal[index + offset + 2] = coef * (temp[i] - 2.0 * OneBySqrt42 * channel_gain[n][0]);
					}
					else{
						Detected_Signal[index + offset + 2] = coef * (-temp[i] - 2.0 * OneBySqrt42 * channel_gain[n][0]);
					}
				}
			}
		}
		index += MODULATION;

	}
	deinterleaver4SD(Detected_Signal, diLLR[data_number - Np]);
}

void CoherentDetection4HD(
	int data_number,
	double(*DFT_signal)[N][2],
	double(*W)[BS][CL][2],
	int(*diLLR)[Nd][N_cbps],
	int(*Decision)){

	int ts, n, index, Detected_Signal[KI][N_cbps];
	double e[2][2], temp1[2];
	int s3_index, s1_index, m1, m3;
	double s_qpsk[4][2] = { { 0.707106781, OneBySqrt2 }, { OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, OneBySqrt2 } };
	double temp[2];

	index = 0;

	for (n = 0; n < N; n++) {
		temp[0] = temp[1] = 0.0;
		temp1[0] = temp1[1] = 0.0;
		e[0][0] = 0.0; e[0][1] = 0.0;

#ifdef  MMSE	/* MMSE and MRC use different coefficients in their respective functions, therefore here they need algorithm specific modifications. */
		for (ts = 0; ts < TS; ts++) {
			temp[0] += W[n][0][ts][0] * DFT_signal[ts][n][0] - W[n][0][ts][1] * DFT_signal[ts][n][1];
			temp[1] += W[n][0][ts][0] * DFT_signal[ts][n][1] + W[n][0][ts][1] * DFT_signal[ts][n][0];
		}
#endif
#ifdef  RLS	
		for (ts = 0; ts < TS; ts++) {
			temp[0] += W[n][0][ts][0] * DFT_signal[ts][n][0] - W[n][0][ts][1] * DFT_signal[ts][n][1];
			temp[1] += W[n][0][ts][0] * DFT_signal[ts][n][1] + W[n][0][ts][1] * DFT_signal[ts][n][0];
		}
		//		temp[0] = temp[0] / channel_power[index][0];
		//		temp[1] = temp[1] / channel_power[index][0];
#endif
#ifdef  MRC 
		for (ts = 0; ts < TS; ts++) {	/* Recall that for MRC, optimal combining is the channel itself. :: her = h^H */
			temp[0] += her[n][ts][0] * DFT_signal[ts][n][0] - her[n][ts][1] * DFT_signal[ts][n][1];
			temp[1] += her[n][ts][0] * DFT_signal[ts][n][1] + her[n][ts][1] * DFT_signal[ts][n][0];

		}
		temp[0] = temp[0] / channel_power[index][0];
		temp[1] = temp[1] / channel_power[index][0];
#endif

		if (MODULATION == QPSK){
			/************************ Hard QPSK decoding *************
			Detected_Signal[index] = (temp[0] >= 0.0) ? 1 : 0;
			Detected_Signal[index + 1] = (temp[1] >= 0.0) ? 1 : 0; */
			//s1_index = 0; s3_index = 0;
			for (m1 = 0; m1 < 4; m1++){
				for (m3 = 0; m3 < 4; m3++){
					temp1[0] = temp[0] + (W[n][0][CL - 1][0] * s_qpsk[m3][0] - W[n][0][CL - 1][1] * s_qpsk[m3][1]);
					temp1[1] = temp[1] + (W[n][0][CL - 1][0] * s_qpsk[m3][1] + W[n][0][CL - 1][1] * s_qpsk[m3][0]);
					//cout << temp[0] << "\t" << temp[1] << "\n";
					//cout << temp1[0] << "\t" << temp1[1] << "\n";
					//cout << s_qpsk[m3][0] << "\t" << s_qpsk[m3][1] << endl;
					e[1][0] = s_qpsk[m1][0] - temp1[0];
					e[1][1] = s_qpsk[m1][1] - temp1[1];
					//cout << m1 << m3 << "\n";
					//cout << e[1][0] << "\t" << e[1][1] << "\n";
					//_getch();
					if (m1 == 0 && m3 == 0){
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
			//cout << e[0][0] << "\t" << e[0][1] << "\n";
			//cout << s1_index << s3_index << endl;
			//_getch();
			Detected_Signal[0][index] = (s_qpsk[s1_index][0] >= 0.0) ? 1 : 0;
			Detected_Signal[0][index + 1] = (s_qpsk[s1_index][1] >= 0.0) ? 1 : 0;
			Detected_Signal[1][index] = (s_qpsk[s3_index][0] >= 0.0) ? 1 : 0;
			Detected_Signal[1][index + 1] = (s_qpsk[s3_index][1] >= 0.0) ? 1 : 0;
			Decision[0] = s1_index; Decision[1] = s3_index;
			/************************ Hard QAM16 decoding *************/
		}
		else if (MODULATION == QAM16){
			Detected_Signal[0][index] = (temp[0] >= 0.0) ? 1 : 0;
			Detected_Signal[0][index + 1] = (-2.0*OneBySqrt10 <= temp[0] && temp[0] <= 2.0*OneBySqrt10) ? 1 : 0;
			Detected_Signal[0][index + 2] = (temp[1] >= 0.0) ? 1 : 0;
			Detected_Signal[0][index + 3] = (-2.0*OneBySqrt10 <= temp[1] && temp[1] <= 2.0*OneBySqrt10) ? 1 : 0;
		}

		/************************ Hard QAM64 Mod decoding *************/

		else if (MODULATION == QAM64){
			if (temp[0] > 0.0)Detected_Signal[0][index] = 1;
			else Detected_Signal[0][index] = 0;


			/*.....................*/
			if (fabs(temp[0]) > 4.0 / sqrt(42.0)){
				Detected_Signal[0][index + 1] = 0;
				if (fabs(temp[0]) > 6.0 / sqrt(42.0)) Detected_Signal[0][index + 2] = 0;
				else Detected_Signal[0][index + 2] = 1;
			}
			/*.....................*/
			else{
				Detected_Signal[0][index + 1] = 1;
				if (fabs(temp[0]) > 2.0 / sqrt(42.0)) Detected_Signal[0][index + 2] = 1;
				else Detected_Signal[0][index + 2] = 0;
			}

			/************************ ********************** *************/

			if (temp[1] > 0.0) Detected_Signal[0][index + 3] = 1;
			else Detected_Signal[0][index + 3] = 0;

			/*.....................*/
			if (fabs(temp[1]) > 4.0 / sqrt(42.0)){
				Detected_Signal[0][index + 4] = 0;
				if (fabs(temp[1]) > 6.0 / sqrt(42.0)) Detected_Signal[0][index + 5] = 0;
				else Detected_Signal[0][index + 5] = 1;

			}
			/*.....................*/
			else{
				Detected_Signal[0][index + 4] = 1;
				if (fabs(temp[1]) > 2.0 / sqrt(42.0)) Detected_Signal[0][index + 5] = 1;
				else Detected_Signal[0][index + 5] = 0;
			}
		}
		index += MODULATION;
	}
	deinterleaver4HD(Detected_Signal[0], diLLR[0][data_number - Np]);
	deinterleaver4HD(Detected_Signal[1], diLLR[1][data_number - Np]);
}
