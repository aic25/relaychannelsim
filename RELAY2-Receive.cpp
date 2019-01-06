#include	"const.h"
extern double h_r2[PATH][2], CNR;


void Relay2_receiver_af(double(*received_signal)[2], double(*transmitted_signal)[2]){
	int	i, k, time;
	double frequency_signal[Nd][Nf][2], data_signal[SAMPLEN][2], time_signal[Nd][Nf][2], H_r2[N][2], DFT_signal[PACKETN][N][2], temp;
	double beta2[PACKETN];
	temp = 0;

#ifdef CUNNING_Relay
	Relay_ideal_channel_estimation(H_r2);	/*  */
#else
#endif
	/*for (size_t i = 0; i < N; i++)
	{
	temp += (sqr(H_r2[i][0]) + sqr(H_r2[i][0])) / N;
	}
	beta2[0] = 1 / sqrt(temp);*/

	for (i = Np; i < PACKETN; i++){
		Relay_DFT(i, received_signal, DFT_signal[i]);	/* Switch to frequency domain considering the order in total data sent. */

		for (size_t j = 0; j < N; j++)
		{
			temp += (sqr(DFT_signal[i][j][0]) + sqr(DFT_signal[i][j][1])) / N;
		}
		beta2[i] = 1 / sqrt(temp);
		temp = 0;

		for (size_t j = 0; j < N; j++)
		{
			DFT_signal[i][j][0] = beta2[i] * DFT_signal[i][j][0];
			DFT_signal[i][j][1] = beta2[i] * DFT_signal[i][j][1];
		}
		for (size_t n = 0; n < Nf; n++) frequency_signal[i - Np][n][0] = frequency_signal[i - Np][n][1] = 0.0;

		for (size_t n = 0; n < N; n++) {
			frequency_signal[i - Np][(n + Nf - N / 2) % Nf][0] = DFT_signal[i][n][0];
			frequency_signal[i - Np][(n + Nf - N / 2) % Nf][1] = DFT_signal[i][n][1];
		}
	}
	for (size_t nd = 0; nd < Nd; nd++){
		IFFT(frequency_signal[nd], time_signal[nd]);
	}
	for (size_t nd = 0; nd < Nd; nd++){
		data_GI_insert(time_signal[nd], data_signal);
		for (size_t t_index = 0; t_index < SAMPLEN; t_index++){
			time = (Np + nd) * SAMPLEN + t_index;
			transmitted_signal[time][0] = data_signal[t_index][0];
			transmitted_signal[time][1] = data_signal[t_index][1];
		}
	}
}

void Relay2_receiver(double(*received_signal)[2], int(*rbit)[N_dbps]){

	int	i, k;
	double H_r2[N][2], DFT_signal[PACKETN][N][2];
	int decode[Nd][N_cbps];

#ifdef RLY_PROPOGATION
#ifdef CUNNING_Relay
	Relay2_ideal_channel_estimation(H_r2);
#else
#endif  
#else
	for (size_t j = 0; j < N; j++)
	{
		H_r2[j][0] = 1;
		H_r2[j][1] = 0;
	}
#endif // RLY_PROPOGATION


	for (i = Np; i < PACKETN; i++){
		Relay_DFT(i, received_signal, DFT_signal[i]);
		Relay_CoherentDetection(i, DFT_signal[i], H_r2, decode);
	}

	for (k = 0; k < Nd; k++){
		HD_Viterbi_decoder(decode[k], rbit[k]);
	}
}

void Relay2_ideal_channel_estimation(double(*H)[2]){
	int  n;
	double Ht[Nf][2], ht[Nf][2], coef;

	coef = sqrt((double)Nf);

	for (n = 0; n < Nf; n++){
		if ((n % (DELAY) == 0) && (n / DELAY) < PATH){
			ht[n][0] = h_r2[n / DELAY][0];
			ht[n][1] = h_r2[n / DELAY][1];
		}
		else{
			ht[n][0] = ht[n][1] = 0.0;
		}
	}
	FFT(ht, Ht);

	for (n = 0; n < N; n++){
		H[n][0] = coef * Ht[(n + Nf - N / 2) % Nf][0];
		H[n][1] = coef * Ht[(n + Nf - N / 2) % Nf][1];
	}
}




