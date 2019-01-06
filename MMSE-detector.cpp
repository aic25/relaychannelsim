#include "const-function-declerations.h"	
extern double CNR, channel_gain[N][2], p_in[N], her[N][TS][2], channel_power[N][2], H[N][TS][BS][2], EigSort[4];
extern int loop, R_index[N], icerik[4];
extern Matrix4cf XCorr[N];
void CoherentDetection4SD(int data_number, double(*DFT_signal)[N][2], double(*WH)[BS][TS][2], double(*diLLR)[N_cbps]){
	int ts, f_index, index, i, offset;
	double temp[2], temp2[2], detected_signal[N_cbps], coef,
		no = pow(10.0, -CNR / 10.0) / CODING_RATE;

	index = 0;
	for (f_index = 0; f_index < N; f_index++) {
		temp[0] = temp[1] = 0.0;
		temp2[0] = temp2[1] = 0.0;

#ifndef MRC
		for (ts = 0; ts < TS; ts++) {

			temp[0] += WH[f_index][0][ts][0] * DFT_signal[ts][f_index][0] - WH[f_index][0][ts][1] * DFT_signal[ts][f_index][1];
			temp[1] += WH[f_index][0][ts][0] * DFT_signal[ts][f_index][1] + WH[f_index][0][ts][1] * DFT_signal[ts][f_index][0];
		}
#else 
		for (ts = 0; ts < TS; ts++) {	/* Recall that for MRC, optimal combining is the channel itself. :: her = h^H */
			temp[0] += her[f_index][ts][0] * DFT_signal[ts][f_index][0] - her[f_index][ts][1] * DFT_signal[ts][f_index][1];
			temp[1] += her[f_index][ts][0] * DFT_signal[ts][f_index][1] + her[f_index][ts][1] * DFT_signal[ts][f_index][0];
		}
#endif  		

		if (MODULATION == QPSK){
			/************************ QPSK ****************************/
			/*
			temp2[0] = (channel_gain[f_index][0] * temp[0] - channel_gain[f_index][1] * temp[1]) / p_in[f_index];
			temp2[1] = (channel_gain[f_index][0] * temp[1] + channel_gain[f_index][1] * temp[0]) / p_in[f_index];

			detected_signal[index] = temp2[0];
			detected_signal[index + 1] = temp2[1];
			*/
			detected_signal[index] = temp[0];
			detected_signal[index + 1] = temp[1];

			/***********************************************************/
		}
		else if (MODULATION == QAM16){
			coef = 2.0*OneBySqrt10 / (1.0 - channel_gain[f_index][0]);
			if (fabs(temp[0]) <= 2.0 * OneBySqrt10 * channel_gain[f_index][0]){
				detected_signal[index] = coef * OneBySqrt10 * temp[0];
			}
			else{
				if (temp[0] >= 0.0){
					detected_signal[index] = coef * OneBySqrt10 * (2.0 * temp[0] - 2.0 * OneBySqrt10 * channel_gain[f_index][0]);
				}
				else {
					detected_signal[index] = coef * OneBySqrt10 * (2.0 * temp[0] + 2.0 * OneBySqrt10 * channel_gain[f_index][0]);
				}
			}
			detected_signal[index + 1] = coef * OneBySqrt10 * (2.0 * OneBySqrt10 * channel_gain[f_index][0] - fabs(temp[0]));

			if (fabs(temp[1]) <= 2.0 * OneBySqrt10 * channel_gain[f_index][0]){
				detected_signal[index + 2] = coef * OneBySqrt10 * temp[1];
			}
			else{
				if (temp[1] >= 0.0){
					detected_signal[index + 2] = coef * OneBySqrt10 * (2.0 * temp[1] - 2.0 * OneBySqrt10 * channel_gain[f_index][0]);
				}
				else {
					detected_signal[index + 2] = coef * OneBySqrt10 * (2.0 * temp[1] + 2.0 * OneBySqrt10 * channel_gain[f_index][0]);
				}
			}
			detected_signal[index + 3] = coef * OneBySqrt10 * (2.0 * OneBySqrt10 * channel_gain[f_index][0] - fabs(temp[1]));


		}
		else if (MODULATION == QAM64){
			coef = 2.0*OneBySqrt42 / (1.0 - channel_gain[f_index][0]);
			for (i = 0; i<2; i++){
				offset = i * 3;
				if (fabs(temp[i]) > 6.0 * OneBySqrt42 * channel_gain[f_index][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset] = coef * (4.0 * temp[i] - 12.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
					else{
						detected_signal[index + offset] = coef * (4.0 * temp[i] + 12.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
				}
				else if (fabs(temp[i]) > 4.0 * OneBySqrt42 * channel_gain[f_index][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset] = coef * (3.0 * temp[i] - 6.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
					else{
						detected_signal[index + offset] = coef * (3.0 * temp[i] + 6.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
				}
				else if (fabs(temp[i]) > 2.0 * OneBySqrt42 * channel_gain[f_index][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset] = coef * (2.0 * temp[i] - 2.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
					else{
						detected_signal[index + offset] = coef * (2.0 * temp[i] + 2.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
				}
				else{
					detected_signal[index + offset] = coef * temp[i];
				}

				if (fabs(temp[i]) > 6.0 * OneBySqrt42 * channel_gain[f_index][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 1] = coef * (-2.0 * temp[i] + 10.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
					else{
						detected_signal[index + offset + 1] = coef * (2.0 * temp[i] + 10.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
				}
				else if (fabs(temp[i]) > 2.0 * OneBySqrt42 * channel_gain[f_index][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 1] = coef * (-temp[i] + 4.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
					else {
						detected_signal[index + offset + 1] = coef * (temp[i] + 4.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
				}
				else{
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 1] = coef * (-2.0 * temp[i] + 6.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
					else{
						detected_signal[index + offset + 1] = coef * (2.0 * temp[i] + 6.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
				}

				if (fabs(temp[i]) > 4.0 * OneBySqrt42 * channel_gain[f_index][0]){
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 2] = coef * (-temp[i] + 6.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
					else{
						detected_signal[index + offset + 2] = coef * (temp[i] + 6.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
				}
				else{
					if (temp[i] >= 0.0){
						detected_signal[index + offset + 2] = coef * (temp[i] - 2.0 * OneBySqrt42 * channel_gain[f_index][0]);
					}
					else{
						detected_signal[index + offset + 2] = coef * (-temp[i] - 2.0 * OneBySqrt42 * channel_gain[f_index][0]);
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
	double(*W)[4][4][2],
	int(*diLLR)[Nd][N_cbps],
	int(*Decision)){

	int ts, f_index, index, detected_signal[KI][N_cbps], s3_index, s1_index, m1, m3;
	double
		index_difference,
		temp[2],
		temp1[2],
		temp2[2],
		temp3[2],
		e[2][2],
		e_1[2],
		e_2[2],
		e_3[2],
		s_qpsk[4][2] = {
			{ OneBySqrt2, OneBySqrt2 },
			{ OneBySqrt2, -OneBySqrt2 },
			{ -OneBySqrt2, -OneBySqrt2 },
			{ -OneBySqrt2, OneBySqrt2 }
	};

	Matrix4cf XXH;
	Matrix4cf Temp4cfM;
	Vector4cf X;
	ComplexEigenSolver<Matrix4cf> Eigs;
	std::ptrdiff_t a;
	//std::complex<float> minEigVal;
	float minEigVal_[2];

	index = 0;

	for (f_index = 0; f_index < N; f_index++) {
		temp[0] = temp[1] = 0.0;
		temp1[0] = temp1[1] = 0.0;
		temp2[0] = temp2[1] = 0.0;
		temp3[0] = temp3[1] = 0.0;
		e[0][0] = 0.0; e[0][1] = 0.0;
		e_1[0] = 0.0; e_1[1] = 0.0;
		e_2[0] = 0.0; e_2[1] = 0.0;
		e_3[0] = 0.0; e_3[1] = 0.0;
		index_difference = data_number - (Np - PREAMBLE_LENGTH) + 1;

#ifdef  MMSE	/* MMSE and MRC use different coefficients in their respective functions, therefore here they need algorithm specific modifications. */
		for (ts = 0; ts < TS; ts++) {
			temp[0] += W[f_index][0][ts][0] * DFT_signal[ts][f_index][0] - W[f_index][0][ts][1] * DFT_signal[ts][f_index][1];
			temp[1] += W[f_index][0][ts][0] * DFT_signal[ts][f_index][1] + W[f_index][0][ts][1] * DFT_signal[ts][f_index][0];
		}
#endif
#ifdef  RLS		 /* not that this is actually necessary */
		for (ts = 0; ts < TS; ts++) {
			temp[0] += W[f_index][0][ts][0] * DFT_signal[ts][f_index][0] + W[f_index][0][ts][1] * DFT_signal[ts][f_index][1];
			temp[1] += W[f_index][0][ts][0] * DFT_signal[ts][f_index][1] - W[f_index][0][ts][1] * DFT_signal[ts][f_index][0];
			if (Metric_Order>1)
			{
				temp1[0] += +W[f_index][(R_index[f_index] + 1) % 2][ts][0] * DFT_signal[ts][f_index][0] + W[f_index][(R_index[f_index] + 1) % 2][ts][1] * DFT_signal[ts][f_index][1];
				temp1[1] += +W[f_index][(R_index[f_index] + 1) % 2][ts][0] * DFT_signal[ts][f_index][1] - W[f_index][(R_index[f_index] + 1) % 2][ts][1] * DFT_signal[ts][f_index][0];
				if (Metric_Order > 2)
				{
					temp2[0] += +W[f_index][2][ts][0] * DFT_signal[ts][f_index][0] + W[f_index][2][ts][1] * DFT_signal[ts][f_index][1];
					temp2[1] += +W[f_index][2][ts][0] * DFT_signal[ts][f_index][1] - W[f_index][2][ts][1] * DFT_signal[ts][f_index][0];
					if (Metric_Order > 3)
					{
						temp3[0] += +W[f_index][3][ts][0] * DFT_signal[ts][f_index][0] + W[f_index][3][ts][1] * DFT_signal[ts][f_index][1];
						temp3[1] += +W[f_index][3][ts][0] * DFT_signal[ts][f_index][1] - W[f_index][3][ts][1] * DFT_signal[ts][f_index][0];
					}
				}
			}
		}
#endif
#ifdef  MRC 
		for (ts = 0; ts < TS; ts++) {	/* Recall that for MRC, optimal combining is the channel itself. :: her = h^H */
			temp[0] += her[f_index][ts][0] * DFT_signal[ts][f_index][0] - her[f_index][ts][1] * DFT_signal[ts][f_index][1];
			temp[1] += her[f_index][ts][0] * DFT_signal[ts][f_index][1] + her[f_index][ts][1] * DFT_signal[ts][f_index][0];
		}
		temp[0] = temp[0] / channel_power[index][0];
		temp[1] = temp[1] / channel_power[index][0];
#endif


		if (MODULATION == QPSK){
			/************************ Hard QPSK decoding *************/
			for (m1 = 0; m1 < 4; m1++){
				for (m3 = 0; m3 < 4; m3++){
					e[1][0] = temp[0] + (+W[f_index][R_index[f_index]][2][0] * s_qpsk[m1][0] + W[f_index][R_index[f_index]][2][1] * s_qpsk[m1][1]) + (+W[f_index][R_index[f_index]][3][0] * s_qpsk[m3][0] + W[f_index][R_index[f_index]][3][1] * s_qpsk[m3][1]);
					e[1][1] = temp[1] + (+W[f_index][R_index[f_index]][2][0] * s_qpsk[m1][1] - W[f_index][R_index[f_index]][2][1] * s_qpsk[m1][0]) + (+W[f_index][R_index[f_index]][3][0] * s_qpsk[m3][1] - W[f_index][R_index[f_index]][3][1] * s_qpsk[m3][0]);
					e[1][0] = (sqr(e[1][0]) + sqr(e[1][1])) / EigSort[icerik[R_index[f_index]]];
					if (Metric_Order > 1)	{
						e_1[0] = temp1[0] + (+W[f_index][(R_index[f_index] + 1) % 2][2][0] * s_qpsk[m1][0] + W[f_index][(R_index[f_index] + 1) % 2][2][1] * s_qpsk[m1][1]) + (+W[f_index][(R_index[f_index] + 1) % 2][3][0] * s_qpsk[m3][0] + W[f_index][(R_index[f_index] + 1) % 2][3][1] * s_qpsk[m3][1]);
						e_1[1] = temp1[1] + (+W[f_index][(R_index[f_index] + 1) % 2][2][0] * s_qpsk[m1][1] - W[f_index][(R_index[f_index] + 1) % 2][2][1] * s_qpsk[m1][0]) + (+W[f_index][(R_index[f_index] + 1) % 2][3][0] * s_qpsk[m3][1] - W[f_index][(R_index[f_index] + 1) % 2][3][1] * s_qpsk[m3][0]);
						e_1[0] = (sqr(e_1[0]) + sqr(e_1[1])) / EigSort[icerik[(R_index[f_index] + 1) % 2]];
						e[1][0] += e_1[0];
						if (Metric_Order > 2){
							e_2[0] = temp2[0] + (+W[f_index][2][2][0] * s_qpsk[m1][0] + W[f_index][2][2][1] * s_qpsk[m1][1]) + (+W[f_index][2][3][0] * s_qpsk[m3][0] + W[f_index][2][3][1] * s_qpsk[m3][1]);
							e_2[1] = temp2[1] + (+W[f_index][2][2][0] * s_qpsk[m1][1] - W[f_index][2][2][1] * s_qpsk[m1][0]) + (+W[f_index][2][3][0] * s_qpsk[m3][1] - W[f_index][2][3][1] * s_qpsk[m3][0]);
							e_2[0] = (sqr(e_2[0]) + sqr(e_2[1])) / EigSort[icerik[2]];
							e[1][0] += e_2[0];
							if (Metric_Order > 3){
								e_3[0] = temp3[0] + (+W[f_index][3][2][0] * s_qpsk[m1][0] + W[f_index][3][2][1] * s_qpsk[m1][1]) + (+W[f_index][3][3][0] * s_qpsk[m3][0] + W[f_index][3][3][1] * s_qpsk[m3][1]);
								e_3[1] = temp3[1] + (+W[f_index][3][2][0] * s_qpsk[m1][1] - W[f_index][3][2][1] * s_qpsk[m1][0]) + (+W[f_index][3][3][0] * s_qpsk[m3][1] - W[f_index][3][3][1] * s_qpsk[m3][0]);
								e_3[0] = (sqr(e_3[0]) + sqr(e_3[1])) / EigSort[icerik[3]];
								e[1][0] += e_3[0];
							}
						}
					}
					if (m1 == 0 && m3 == 0){
						e[0][0] = e[1][0];
						s1_index = 0;
						s3_index = 0;
					}
					else if (sqrt(sqr(e[1][0])) < sqrt(sqr(e[0][0]))){
						e[0][0] = e[1][0];
						s1_index = m1;
						s3_index = m3;
					}
				}
			}
#ifdef EVDF
			X[0] = complex<double>(DFT_signal[0][f_index][0], DFT_signal[0][f_index][1]);
			X[1] = complex<double>(DFT_signal[1][f_index][0], DFT_signal[1][f_index][1]);
			X[2] = complex<double>(s_qpsk[s1_index][0], s_qpsk[s1_index][1]);
			X[3] = complex<double>(s_qpsk[s3_index][0], s_qpsk[s3_index][1]);
			XXH = X*X.adjoint();
			Temp4cfM = XCorr[f_index];
			XCorr[f_index] = (1 / index_difference)*XXH + ((index_difference - 1) / index_difference)*Temp4cfM;
			Eigs.compute(XCorr[f_index]);

			EigSort[0] = Eigs.eigenvalues().real()[0];
			EigSort[1] = Eigs.eigenvalues().real()[1];
			EigSort[2] = Eigs.eigenvalues().real()[2];
			EigSort[3] = Eigs.eigenvalues().real()[3];
			sirala(EigSort, 4, icerik);
			EigSort[0] = Eigs.eigenvalues().real()[0];
			EigSort[1] = Eigs.eigenvalues().real()[1];
			EigSort[2] = Eigs.eigenvalues().real()[2];
			EigSort[3] = Eigs.eigenvalues().real()[3];

			W[f_index][R_index[f_index]][0][0] = Eigs.eigenvectors().col(icerik[0])[0].real();
			W[f_index][R_index[f_index]][0][1] = Eigs.eigenvectors().col(icerik[0])[0].imag();
			W[f_index][R_index[f_index]][1][0] = Eigs.eigenvectors().col(icerik[0])[1].real();
			W[f_index][R_index[f_index]][1][1] = Eigs.eigenvectors().col(icerik[0])[1].imag();
			W[f_index][R_index[f_index]][2][0] = Eigs.eigenvectors().col(icerik[0])[2].real();
			W[f_index][R_index[f_index]][2][1] = Eigs.eigenvectors().col(icerik[0])[2].imag();
			W[f_index][R_index[f_index]][3][0] = Eigs.eigenvectors().col(icerik[0])[3].real();
			W[f_index][R_index[f_index]][3][1] = Eigs.eigenvectors().col(icerik[0])[3].imag();

			W[f_index][(R_index[f_index] + 1) % 2][0][0] = Eigs.eigenvectors().col(icerik[1])[0].real();
			W[f_index][(R_index[f_index] + 1) % 2][0][1] = Eigs.eigenvectors().col(icerik[1])[0].imag();
			W[f_index][(R_index[f_index] + 1) % 2][1][0] = Eigs.eigenvectors().col(icerik[1])[1].real();
			W[f_index][(R_index[f_index] + 1) % 2][1][1] = Eigs.eigenvectors().col(icerik[1])[1].imag();
			W[f_index][(R_index[f_index] + 1) % 2][2][0] = Eigs.eigenvectors().col(icerik[1])[2].real();
			W[f_index][(R_index[f_index] + 1) % 2][2][1] = Eigs.eigenvectors().col(icerik[1])[2].imag();
			W[f_index][(R_index[f_index] + 1) % 2][3][0] = Eigs.eigenvectors().col(icerik[1])[3].real();
			W[f_index][(R_index[f_index] + 1) % 2][3][1] = Eigs.eigenvectors().col(icerik[1])[3].imag();

			W[f_index][2][0][0] = Eigs.eigenvectors().col(icerik[2])[0].real();
			W[f_index][2][0][1] = Eigs.eigenvectors().col(icerik[2])[0].imag();
			W[f_index][2][1][0] = Eigs.eigenvectors().col(icerik[2])[1].real();
			W[f_index][2][1][1] = Eigs.eigenvectors().col(icerik[2])[1].imag();
			W[f_index][2][2][0] = Eigs.eigenvectors().col(icerik[2])[2].real();
			W[f_index][2][2][1] = Eigs.eigenvectors().col(icerik[2])[2].imag();
			W[f_index][2][3][0] = Eigs.eigenvectors().col(icerik[2])[3].real();
			W[f_index][2][3][1] = Eigs.eigenvectors().col(icerik[2])[3].imag();

			W[f_index][3][0][0] = Eigs.eigenvectors().col(icerik[3])[0].real();
			W[f_index][3][0][1] = Eigs.eigenvectors().col(icerik[3])[0].imag();
			W[f_index][3][1][0] = Eigs.eigenvectors().col(icerik[3])[1].real();
			W[f_index][3][1][1] = Eigs.eigenvectors().col(icerik[3])[1].imag();
			W[f_index][3][2][0] = Eigs.eigenvectors().col(icerik[3])[2].real();
			W[f_index][3][2][1] = Eigs.eigenvectors().col(icerik[3])[2].imag();
			W[f_index][3][3][0] = Eigs.eigenvectors().col(icerik[3])[3].real();
			W[f_index][3][3][1] = Eigs.eigenvectors().col(icerik[3])[3].imag();
#endif
			detected_signal[0][index] = (R_index[f_index] == 0) ? ((s_qpsk[s1_index][0] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][0] >= 0.0) ? 1 : 0);
			detected_signal[0][index + 1] = (R_index[f_index] == 0) ? ((s_qpsk[s1_index][1] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][1] >= 0.0) ? 1 : 0);
			detected_signal[1][index] = (R_index[f_index] == 1) ? ((s_qpsk[s1_index][0] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][0] >= 0.0) ? 1 : 0);
			detected_signal[1][index + 1] = (R_index[f_index] == 1) ? ((s_qpsk[s1_index][1] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][1] >= 0.0) ? 1 : 0);
			Decision[0] = (R_index[f_index] == 0) ? s1_index : s3_index; Decision[1] = (R_index[f_index] == 1) ? s1_index : s3_index;
			/************************ Hard QAM16 decoding *************/
		}
		else if (MODULATION == QAM16){
			detected_signal[R_index[f_index]][index] = (temp[0] >= 0.0) ? 1 : 0;
			detected_signal[R_index[f_index]][index + 1] = (-2.0*OneBySqrt10 <= temp[0] && temp[0] <= 2.0*OneBySqrt10) ? 1 : 0;
			detected_signal[R_index[f_index]][index + 2] = (temp[1] >= 0.0) ? 1 : 0;
			detected_signal[R_index[f_index]][index + 3] = (-2.0*OneBySqrt10 <= temp[1] && temp[1] <= 2.0*OneBySqrt10) ? 1 : 0;
		}

		/************************ Hard QAM64 Mod decoding *************/

		else if (MODULATION == QAM64){
			if (temp[0] > 0.0)detected_signal[R_index[f_index]][index] = 1;
			else detected_signal[R_index[f_index]][index] = 0;


			/*.....................*/
			if (fabs(temp[0]) > 4.0 / sqrt(42.0)){
				detected_signal[R_index[f_index]][index + 1] = 0;
				if (fabs(temp[0]) > 6.0 / sqrt(42.0)) detected_signal[R_index[f_index]][index + 2] = 0;
				else detected_signal[R_index[f_index]][index + 2] = 1;
			}
			/*.....................*/
			else{
				detected_signal[R_index[f_index]][index + 1] = 1;
				if (fabs(temp[0]) > 2.0 / sqrt(42.0)) detected_signal[R_index[f_index]][index + 2] = 1;
				else detected_signal[R_index[f_index]][index + 2] = 0;
			}

			/************************ ********************** *************/

			if (temp[1] > 0.0) detected_signal[R_index[f_index]][index + 3] = 1;
			else detected_signal[R_index[f_index]][index + 3] = 0;

			/*.....................*/
			if (fabs(temp[1]) > 4.0 / sqrt(42.0)){
				detected_signal[R_index[f_index]][index + 4] = 0;
				if (fabs(temp[1]) > 6.0 / sqrt(42.0)) detected_signal[R_index[f_index]][index + 5] = 0;
				else detected_signal[R_index[f_index]][index + 5] = 1;

			}
			/*.....................*/
			else{
				detected_signal[R_index[f_index]][index + 4] = 1;
				if (fabs(temp[1]) > 2.0 / sqrt(42.0)) detected_signal[R_index[f_index]][index + 5] = 1;
				else detected_signal[R_index[f_index]][index + 5] = 0;
			}
		}
		index += MODULATION;
	}
	deinterleaver4HD(detected_signal[0], diLLR[0][data_number - Np]);
	deinterleaver4HD(detected_signal[1], diLLR[1][data_number - Np]);
}
