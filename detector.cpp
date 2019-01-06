#include "const.h"
extern double channel_gain[N][2], p_in[N], her[N][TS][2], channel_power[N][2], EigSort[4];
extern int loop, icerik[4];
extern Matrix4cf XCorr[N];
//extern float minEigVal;

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
	int(*R_index),
	int data_number,
	double(*DFT_signal)[N][2],
	double(*W)[4][4][2],
	int(*diLLR)[Nd][N_cbps],
	int(*Decision)){

	int ts, n, index, Detected_Signal[KI][N_cbps];
	double e[2][2], temp1[2], e_1[2], e_2[2], e_3[2];
	int s3_index, s1_index, m1, m3, a1, a3;
	double s_qpsk[4][2] = { { 0.707106781, OneBySqrt2 }, { OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, OneBySqrt2 } };
	double temp[2], temp_1[2], temp_2[2], temp_3[2];

	double index_difference;
	Matrix4cf XXH;
	Vector4cf X;
	ComplexEigenSolver<Matrix4cf> Eigs;
	std::ptrdiff_t a;
	//std::complex<float> minEigVal;
	float minEigVal_[2];

	index = 0;

	for (n = 0; n < N; n++) {
		temp[0] = temp[1] = 0.0;
		temp_1[0] = temp_1[1] = 0.0;
		temp_2[0] = temp_2[1] = 0.0;
		temp_3[0] = temp_3[1] = 0.0;
		temp1[0] = temp1[1] = 0.0;
		e[0][0] = 0.0; e[0][1] = 0.0;
		e_1[0] = 0.0; e_1[1] = 0.0;
		e_2[0] = 0.0; e_2[1] = 0.0;
		e_3[0] = 0.0; e_3[1] = 0.0;
		index_difference = data_number - (Np - PREAMBLE_LENGTH) + 1;

#ifdef  MMSE	/* MMSE and MRC use different coefficients in their respective functions, therefore here they need algorithm specific modifications. */
		for (ts = 0; ts < TS; ts++) {
			temp[0] += W[n][0][ts][0] * DFT_signal[ts][n][0] - W[n][0][ts][1] * DFT_signal[ts][n][1];
			temp[1] += W[n][0][ts][0] * DFT_signal[ts][n][1] + W[n][0][ts][1] * DFT_signal[ts][n][0];
		}
#endif
#ifdef  RLS	
		for (ts = 0; ts < TS; ts++) {
			temp[0] += +W[n][R_index[n]][ts][0] * DFT_signal[ts][n][0] + W[n][R_index[n]][ts][1] * DFT_signal[ts][n][1];
			temp[1] += +W[n][R_index[n]][ts][0] * DFT_signal[ts][n][1] - W[n][R_index[n]][ts][1] * DFT_signal[ts][n][0];		  
			if (Metric_Order>1)
			{
				temp_1[0] += +W[n][(R_index[n] + 1) % 2][ts][0] * DFT_signal[ts][n][0] + W[n][(R_index[n] + 1) % 2][ts][1] * DFT_signal[ts][n][1];
				temp_1[1] += +W[n][(R_index[n] + 1) % 2][ts][0] * DFT_signal[ts][n][1] - W[n][(R_index[n] + 1) % 2][ts][1] * DFT_signal[ts][n][0];
				if (Metric_Order>2)
				{		  
					temp_2[0] += +W[n][2][ts][0] * DFT_signal[ts][n][0] + W[n][2][ts][1] * DFT_signal[ts][n][1];
					temp_2[1] += +W[n][2][ts][0] * DFT_signal[ts][n][1] - W[n][2][ts][1] * DFT_signal[ts][n][0];		 
					if (Metric_Order>3)
					{
						temp_3[0] += +W[n][3][ts][0] * DFT_signal[ts][n][0] + W[n][3][ts][1] * DFT_signal[ts][n][1];
						temp_3[1] += +W[n][3][ts][0] * DFT_signal[ts][n][1] - W[n][3][ts][1] * DFT_signal[ts][n][0];
					}																								 
				}	
			}																							  
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
			for (m1 = 0; m1 < 4; m1++){
			for (m3 = 0; m3 < 4; m3++){
			X[0].real(DFT_signal[0][n][0]);	X[0].imag(DFT_signal[0][n][1]);
			X[1].real(DFT_signal[1][n][0]);	X[1].imag(DFT_signal[1][n][1]);
			X[2].real(s_qpsk[m1][0]);	X[2].imag(s_qpsk[m1][1]);
			X[3].real(s_qpsk[m3][0]);	X[3].imag(s_qpsk[m3][1]);
			XXH = XCorr + X*X.adjoint();
			Eigs.compute(XXH);
			minEigVal_[1] = Eigs.eigenvalues().real().minCoeff(&a);
			//cout << minEigVal_[0] << " " << minEigVal_[1] << " " << m1 << " " << m3 << endl;	getchar();
			if (m1 == 0 && m3 == 0){
			minEigVal_[0] = minEigVal_[1];
			s1_index = 0;
			s3_index = 0;
			}
			else if (minEigVal_[1] < minEigVal_[0]){
			minEigVal_[0] = minEigVal_[1];
			s1_index = m1;
			s3_index = m3;
			}
			}
			}	   */
			/*	*/
			for (m1 = 0; m1 < 4; m1++){
				for (m3 = 0; m3 < 4; m3++){
					e[1][0] = temp[0] + (+W[n][R_index[n]][2][0] * s_qpsk[m3][0] + W[n][R_index[n]][2][1] * s_qpsk[m3][1]) + (+W[n][R_index[n]][3][0] * s_qpsk[m1][0] + W[n][R_index[n]][3][1] * s_qpsk[m1][1]);
					e[1][1] = temp[1] + (+W[n][R_index[n]][2][0] * s_qpsk[m3][1] - W[n][R_index[n]][2][1] * s_qpsk[m3][0]) + (+W[n][R_index[n]][3][0] * s_qpsk[m1][1] - W[n][R_index[n]][3][1] * s_qpsk[m1][0]);
					e[1][0] = (sqr(e[1][0]) + sqr(e[1][1]))/EigSort[icerik[R_index[n]]];
					//e[1][1] /= EigSort[icerik[R_index[n]]];
					if (Metric_Order > 1)	{
						e_1[0] = temp_1[0] + (+W[n][(R_index[n] + 1) % 2][2][0] * s_qpsk[m3][0] + W[n][(R_index[n] + 1) % 2][2][1] * s_qpsk[m3][1]) + (+W[n][(R_index[n] + 1) % 2][3][0] * s_qpsk[m1][0] + W[n][(R_index[n] + 1) % 2][3][1] * s_qpsk[m1][1]);
						e_1[1] = temp_1[1] + (+W[n][(R_index[n] + 1) % 2][2][0] * s_qpsk[m3][1] - W[n][(R_index[n] + 1) % 2][2][1] * s_qpsk[m3][0]) + (+W[n][(R_index[n] + 1) % 2][3][0] * s_qpsk[m1][1] - W[n][(R_index[n] + 1) % 2][3][1] * s_qpsk[m1][0]);
						e_1[0] = (sqr(e_1[0]) + sqr(e_1[1])) / EigSort[icerik[(R_index[n] + 1) % 2]];
						//e_1[1] /= EigSort[icerik[(R_index[n] + 1) % 2]];
						e[1][0] += e_1[0];
						//e[1][1] += e_1[1];
						if (Metric_Order > 2){
							e_2[0] = temp_2[0] + (+W[n][2][2][0] * s_qpsk[m3][0] + W[n][2][2][1] * s_qpsk[m3][1]) + (+W[n][2][3][0] * s_qpsk[m1][0] + W[n][2][3][1] * s_qpsk[m1][1]);
							e_2[1] = temp_2[1] + (+W[n][2][2][0] * s_qpsk[m3][1] - W[n][2][2][1] * s_qpsk[m3][0]) + (+W[n][2][3][0] * s_qpsk[m1][1] - W[n][2][3][1] * s_qpsk[m1][0]);
							e_2[0] = (sqr(e_2[0]) + sqr(e_2[1])) / EigSort[icerik[2]];
							//e_2[1] /= EigSort[icerik[2]];
							e[1][0] += e_2[0];
							//e[1][1] += e_2[1];
							if (Metric_Order > 3){
								e_3[0] = temp_3[0] + (+W[n][3][2][0] * s_qpsk[m3][0] + W[n][3][2][1] * s_qpsk[m3][1]) + (+W[n][3][3][0] * s_qpsk[m1][0] + W[n][3][3][1] * s_qpsk[m1][1]);
								e_3[1] = temp_3[1] + (+W[n][3][2][0] * s_qpsk[m3][1] - W[n][3][2][1] * s_qpsk[m3][0]) + (+W[n][3][3][0] * s_qpsk[m1][1] - W[n][3][3][1] * s_qpsk[m1][0]);
								e_3[0] = (sqr(e_3[0]) + sqr(e_3[1])) / EigSort[icerik[3]];
								//e_2[1] /= EigSort[icerik[2]];
								e[1][0] += e_3[0];
								//e[1][1] += e_2[1];
							}
						}
					}

					//cout << temp[0] << "\t" << temp[1] << "\n";
					//cout << temp1[0] << "\t" << temp1[1] << "\n";
					//cout << s_qpsk[m3][0] << "\t" << s_qpsk[m3][1] << endl;
					//cout << m1 << m3 << "\n";
					//cout << e[1][0] << "\t" << e[1][1] << "\t" << sqrt(sqr(e[1][0]) + sqr(e[1][1])) << "\n";
					//_getch();
					if (m1 == 0 && m3 == 0){
						e[0][0] = e[1][0];
						//e[0][1] = e[1][1];
						//a1 = 0;
						//a3 = 0;
						s1_index = 0;
						s3_index = 0;
					}
					//else if (sqrt(sqr(e[1][0]) + sqr(e[1][1])) < sqrt(sqr(e[0][0]) + sqr(e[0][1]))){
					else if (sqrt(sqr(e[1][0])) < sqrt(sqr(e[0][0]))){
						e[0][0] = e[1][0];
						e[0][1] = e[1][1];
						//a1 = m1;
						//a3 = m3;
						s1_index = m1;
						s3_index = m3;
					}
				}
			}
			//cout << s1_index << " " << s3_index << endl << a1 << " " << a3 << endl;	getchar();	   

			/*for (m1 = 0; m1 < 4; m1++){
				for (m3 = 0; m3 < 4; m3++){
				e[1][0] = temp[0] + (+W[n][R_index[n]][2][0] * s_qpsk[m3][0] + W[n][R_index[n]][2][1] * s_qpsk[m3][1]) + (+W[n][R_index[n]][3][0] * s_qpsk[m1][0] + W[n][R_index[n]][3][1] * s_qpsk[m1][1]);
				e[1][1] = temp[1] + (+W[n][R_index[n]][2][0] * s_qpsk[m3][1] - W[n][R_index[n]][2][1] * s_qpsk[m3][0]) + (+W[n][R_index[n]][3][0] * s_qpsk[m1][1] - W[n][R_index[n]][3][1] * s_qpsk[m1][0]);
				//cout << temp[0] << "\t" << temp[1] << "\n";
				//cout << temp1[0] << "\t" << temp1[1] << "\n";
				//cout << s_qpsk[m3][0] << "\t" << s_qpsk[m3][1] << endl;
				//cout << m1 << m3 << "\n";
				//cout << e[1][0] << "\t" << e[1][1] << "\t" << sqrt(sqr(e[1][0]) + sqr(e[1][1])) << "\n";
				//_getch();
				if (m1 == 0 && m3 == 0){
				e[0][0] = e[1][0] - minEigVal;
				e[0][1] = e[1][1] - minEigVal;
				//a1 = 0;
				//a3 = 0;
				s1_index = 0;
				s3_index = 0;
				}
				else if (sqrt(sqr(e[1][0] - minEigVal) + sqr(e[1][1] - minEigVal)) < sqrt(sqr(e[0][0]) + sqr(e[0][1]))){
				e[0][0] = e[1][0] - minEigVal;
				e[0][1] = e[1][1] - minEigVal;
				//a1 = m1;
				//a3 = m3;
				s1_index = m1;
				s3_index = m3;
				}
				}
				}	*/

			/*  */
#ifdef DDCU
			/*X[0].real(DFT_signal[0][n][0]);	X[0].imag(DFT_signal[0][n][1]);
			X[1].real(DFT_signal[1][n][0]);	X[1].imag(DFT_signal[1][n][1]);
			X[2].real(s_qpsk[s1_index][0]);	X[2].imag(s_qpsk[s1_index][1]);
			X[3].real(s_qpsk[s3_index][0]);	X[3].imag(s_qpsk[s3_index][1]);		   */
			X[0] = complex<double>(DFT_signal[0][n][0], DFT_signal[0][n][1]);
			X[1] = complex<double>(DFT_signal[1][n][0], DFT_signal[1][n][1]);
			X[2] = complex<double>(s_qpsk[s3_index][0], s_qpsk[s3_index][1]);
			X[3] = complex<double>(s_qpsk[s1_index][0], s_qpsk[s1_index][1]);
			XCorr[n] = (1 / index_difference)*(X*X.adjoint()) + ((index_difference - 1) / index_difference)*XCorr[n];
			/*XCorr[n](2, 2).real(1); XCorr[n](2, 2).imag(0);
			XCorr[n](2, 3).real(0); XCorr[n](2, 2).imag(0);
			XCorr[n](3, 2).real(0); XCorr[n](2, 2).imag(0);
			XCorr[n](3, 3).real(1); XCorr[n](2, 2).imag(0);	*/
			/*XCorr[n](2, 2) = complex<double>(1, 0);
			XCorr[n](2, 3) = complex<double>(0, 0);
			XCorr[n](3, 2) = complex<double>(0, 0);
			XCorr[n](3, 3) = complex<double>(1, 0);		 */
			/*if (n == 5 && data_number==(PACKETN-1)){
				cout << "XCorr: " << endl << XCorr[n] << endl
				<< " XCorr(3,3): " << XCorr[n](2, 2) << endl
				<< " XCorr(3,4): " << XCorr[n](2, 3) << endl
				<< " XCorr(4,3): " << XCorr[n](3, 2) << endl
				<< " XCorr(4,4): " << XCorr[n](3, 3) << endl
				<< endl;
				getchar();
				}		 */
			Eigs.compute(XCorr[n]);
			//minEigVal = Eigs.eigenvalues().real().minCoeff(&a);
			W[n][R_index[n]][0][0] = Eigs.eigenvectors().col(a)[0].real();	W[n][R_index[n]][0][1] = Eigs.eigenvectors().col(a)[0].imag();
			W[n][R_index[n]][1][0] = Eigs.eigenvectors().col(a)[1].real();	W[n][R_index[n]][1][1] = Eigs.eigenvectors().col(a)[1].imag();
			W[n][R_index[n]][2][0] = Eigs.eigenvectors().col(a)[2].real();	W[n][R_index[n]][2][1] = Eigs.eigenvectors().col(a)[2].imag();
			W[n][R_index[n]][3][0] = Eigs.eigenvectors().col(a)[3].real();	W[n][R_index[n]][3][1] = Eigs.eigenvectors().col(a)[3].imag();
#endif // DDCU

			//if (n == 1)cout << minEigVal << " " << a << endl;		  

			//cout << e[0][0] << "\t" << e[0][1] << "\t" << sqrt(sqr(e[0][0]) + sqr(e[0][1])) << "\n";
			//cout << s1_index << s3_index << endl;
			//_getch();
			Detected_Signal[0][index] = (R_index[n] == 0) ? ((s_qpsk[s1_index][0] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][0] >= 0.0) ? 1 : 0);
			Detected_Signal[0][index + 1] = (R_index[n] == 0) ? ((s_qpsk[s1_index][1] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][1] >= 0.0) ? 1 : 0);
			Detected_Signal[1][index] = (R_index[n] == 1) ? ((s_qpsk[s1_index][0] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][0] >= 0.0) ? 1 : 0);
			Detected_Signal[1][index + 1] = (R_index[n] == 1) ? ((s_qpsk[s1_index][1] >= 0.0) ? 1 : 0) : ((s_qpsk[s3_index][1] >= 0.0) ? 1 : 0);
			Decision[0] = (R_index[n] == 0) ? s1_index : s3_index; Decision[1] = (R_index[n] == 1) ? s1_index : s3_index;
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
