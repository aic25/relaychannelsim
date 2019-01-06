#include	"const.h"
extern double h_r[BS][PATH][2];
extern double CNR;

void Relay_receiver(double(*Received_Signal)[2], int(*rbit)[N_dbps], int bs){
	int	i, k;
	double H_r[N][2], DFT_signal[PACKETN][N][2];
	int decode[Nd][N_cbps];

#ifdef RLY_PROPOGATION
#ifdef CUNNING_Relay
	Relay_ideal_channel_estimation(H_r, bs);	/*  */
#else
	/* not that relay need channel estimation? */
#endif 
#else
	for (size_t j = 0; j < N; j++)
	{	 
		H_r[j][0] = 1;
		H_r[j][1] = 0;
	}
#endif // RLY_PROPOGATION


	for (i = Np; i < PACKETN; i++){
		Relay_DFT(i, Received_Signal, DFT_signal[i]);	/* Switch to frequency domain considering the order in total data sent. */
		Relay_CoherentDetection(i, DFT_signal[i], H_r, decode);	/* Relay and receiver uses the same modulation. We can change that if necessary. */
	}

	for (k = 0; k < Nd; k++){
		HD_Viterbi_decoder(decode[k], rbit[k]);	/* ? No soft decision for relay */
	}
}

void Relay_ideal_channel_estimation(double(*H_r)[2], int bs){
	int  n;
	double Ht[Nf][2], ht[Nf][2], coef;

	coef = sqrt((double)Nf);

	for (n = 0; n < Nf; n++){
		if ((n % (DELAY) == 0) && (n / DELAY) < PATH){
			ht[n][0] = h_r[bs][n / DELAY][0];
			ht[n][1] = h_r[bs][n / DELAY][1];
		}
		else{
			ht[n][0] = ht[n][1] = 0.0;
		}
	}
	FFT(ht, Ht);

	for (n = 0; n < N; n++){
		H_r[n][0] = coef * Ht[(n + Nf - N / 2) % Nf][0];
		H_r[n][1] = coef * Ht[(n + Nf - N / 2) % Nf][1];
	}
}

void Relay_DFT(int data_number, double(*Received_Signal)[2], double(*DFT_signal)[2]){
	int	f_index, time, time_index, carrier_index;
	double time_signal[Nf][2], frequency_signal[Nf][2], symbol[N][2];

	for (time_index = 0; time_index < Nf; time_index++) {
		time = data_number * SAMPLEN + GI + time_index;
		time_signal[time_index][0] = Received_Signal[time][0];
		time_signal[time_index][1] = Received_Signal[time][1];
	}

	FFT(time_signal, frequency_signal);

	for (carrier_index = 0; carrier_index < N; carrier_index++) {
		symbol[carrier_index][0] = frequency_signal[(carrier_index + Nf - N / 2) % Nf][0];
		symbol[carrier_index][1] = frequency_signal[(carrier_index + Nf - N / 2) % Nf][1];
	}

	for (f_index = 0; f_index < N; f_index++){
		DFT_signal[f_index][0] = symbol[f_index][0];
		DFT_signal[f_index][1] = symbol[f_index][1];
	}
}

void Relay_CoherentDetection(int data_number, double(*DFT_signal)[2], double(*H_r)[2], int(*decode)[N_cbps]){
	int n, index, Detected_Signal[N_cbps];
	double temp[2], power;
	index = 0;
	for (n = 0; n < N; n++) {
		temp[0] = temp[1] = 0.0;

		temp[0] += H_r[n][0] * DFT_signal[n][0] + H_r[n][1] * DFT_signal[n][1];
		temp[1] += H_r[n][0] * DFT_signal[n][1] - H_r[n][1] * DFT_signal[n][0];

		power = sqr(H_r[n][0]) + sqr(H_r[n][1]);
		temp[0] /= power;
		temp[1] /= power;

		/************************ Hard QPSK decoding *************/
		if (MODULATION == QPSK){
			Detected_Signal[index] = (temp[0] >= 0.0) ? 1 : 0;
			Detected_Signal[index + 1] = (temp[1] >= 0.0) ? 1 : 0;
		}
		/************************ Hard QAM16 decoding *************/
		else if (MODULATION == QAM16){
			Detected_Signal[index] = (temp[0] >= 0.0) ? 1 : 0;
			Detected_Signal[index + 1] = (-2.0*OneBySqrt10 <= temp[0] && temp[0] <= 2.0*OneBySqrt10) ? 1 : 0;
			Detected_Signal[index + 2] = (temp[1] >= 0.0) ? 1 : 0;
			Detected_Signal[index + 3] = (-2.0*OneBySqrt10 <= temp[1] && temp[1] <= 2.0*OneBySqrt10) ? 1 : 0;
		}
		/************************ Hard QAM64 Mod decoding *************/
		else if (MODULATION == QAM64){
			if (temp[0] > 0.0)Detected_Signal[index] = 1;
			else Detected_Signal[index] = 0;

			/*.....................*/
			if (fabs(temp[0]) > 4.0 / sqrt(42.0)){
				Detected_Signal[index + 1] = 0;
				if (fabs(temp[0]) > 6.0 / sqrt(42.0)) Detected_Signal[index + 2] = 0;
				else Detected_Signal[index + 2] = 1;
			}
			/*.....................*/
			else{
				Detected_Signal[index + 1] = 1;
				if (fabs(temp[0]) > 2.0 / sqrt(42.0)) Detected_Signal[index + 2] = 1;
				else Detected_Signal[index + 2] = 0;
			}

			/************************ ********************** *************/

			if (temp[1] > 0.0) Detected_Signal[index + 3] = 1;
			else Detected_Signal[index + 3] = 0;

			/*.....................*/
			if (fabs(temp[1]) > 4.0 / sqrt(42.0)){
				Detected_Signal[index + 4] = 0;
				if (fabs(temp[1]) > 6.0 / sqrt(42.0)) Detected_Signal[index + 5] = 0;
				else Detected_Signal[index + 5] = 1;

			}
			/*.....................*/
			else{
				Detected_Signal[index + 4] = 1;   //bit2
				if (fabs(temp[1]) > 2.0 / sqrt(42.0)) Detected_Signal[index + 5] = 1;  //bit3
				else Detected_Signal[index + 5] = 0;
			}

		}
		index += MODULATION;
	}
	deinterleaver4HD(Detected_Signal, decode[data_number - Np]);
}

