#include	"const.h"
extern double h_r2[PATH][2], CNR;

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




