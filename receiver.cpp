#include "const.h"

extern double transmitted_sigmat[Rx][Tx][BURST][2], mse_sum[LOOPN][NUMBER_OF_STEPS];
extern double wh_sequence[2][Np][N], channel_power[N][2], her[N][TS][2];
double h_dsrd[N][TS][2];
double H[N][Rx][Tx][2];
extern int loop;

void receiver(double(*received_signal)[BURST][2], int(*rbit)[N_dbps]){
	int index, rx, k, f_index;
	int coded_bit[N_cbps];
	double received_vector[FILTER_ORDER + 1][2], DFT_signal[PACKETN][Rx][N][2], e_freqav[2], e[2];
	e_freqav[0] = 0; e_freqav[1] = 0;
	double DFT_signal_ideal[PACKETN][Rx][N][2];	/* Ideal one is for RLS. */

	double
		W_mmse[N][Tx][Rx][2],
		W[N][Tx][Rx][2];

#ifdef SD
	double diLLR[Nd][N_cbps];	/* Pre-estimation buffer? */
#endif
#ifdef HD
	int diLLR[Nd][N_cbps];	/* Pre-estimation buffer? */
#endif 			  
	/************************** CHANNEL ESTIMATION *************************/
#ifdef CUNNING
	ideal_channel_estimation4all_ch(H);
#endif
#ifndef CUNNING
#ifdef ZC
	channel_estimation4ZC_4all_ch(received_signal, H);
#else	/* For the case of WH */
	for (index = 0; index < Np; index++){
		for (rx = 0; rx < Rx; rx++) {
			DFT(rx, index, received_signal, DFT_signal[index]);
		}
		DFT(0, index, transmitted_sigmat[0], DFT_signal_ideal[index]);
	}
	WH_channel_estimation(DFT_signal, H, wh_sequence);
#endif  
#endif // !CUNNING
	h_desired(H, h_dsrd);
	VectorHermite(h_dsrd, her);
	VectorMulVectorToscalar(her, h_dsrd, channel_power);
	/*********************** Weight Generation *****************************/
#ifdef MMSE
	MMSE_Weight_Matrix(H, W);
#endif
#ifdef MRC
	mrc_weight(H);
#endif
#ifdef RLS	/* more to add */
#ifndef Wh
	for (index = 0; index < Np; index++){
		for (rx = 0; rx < Rx; rx++) {
			DFT(rx, index, received_signal, DFT_signal[index]);
		}
		DFT(0, index, transmitted_sigmat[0], DFT_signal_ideal[index]);
	}
#endif // 			
	RLS_initial_weigth(0, DFT_signal, DFT_signal_ideal, W, H);						  
#endif	// RLS
	/****************************** DETECTION ******************************/
#ifdef SD
	double frequency_signal[N][2];
	for (index = Np; index < PACKETN; index++){
		for (rx = 0; rx < Rx; rx++) {
			DFT(rx, index, received_signal, DFT_signal[index]);
		}
		CoherentDetection4SD(index, DFT_signal[index], W, diLLR);
		k = index - Np;
		SD_Viterbi_decoder(diLLR[k], rbit[k]);

#ifdef DDWE
		convolutional_encoder(rbit[k], coded_bit);	/* k is the symbol index */
		interleaver(coded_bit, coded_bit);
		data_insert_trsym(k, coded_bit, frequency_signal);
		for (f_index = 0; f_index < N; f_index++){
			for (rx = 0; rx < Rx; rx++){
				received_vector[rx][0] = DFT_signal[index][rx][f_index][0];
				received_vector[rx][1] = DFT_signal[index][rx][f_index][1];
			}
			DDWE_Update(received_vector, frequency_signal[f_index], W[f_index][0], W[f_index][0], H[f_index], e, f_index);
		}
#endif // DDWE
	}			   
#endif

#ifdef HD
	double frequency_signal[N][2];
	for (index = Np; index < PACKETN; index++){
		for (rx = 0; rx < Rx; rx++) {
			DFT(rx, index, received_signal, DFT_signal[index]);
		}
		CoherentDetection4HD(index, DFT_signal[index], W, diLLR);
		k = index - Np;
		HD_Viterbi_decoder(diLLR[k], rbit[k]);

#ifdef DDWE
		convolutional_encoder(rbit[k], coded_bit);	/* k is the symbol index */
		interleaver(coded_bit, coded_bit);
		data_insert_trsym(k, coded_bit, frequency_signal);
		for (f_index = 0; f_index < N; f_index++){
			for (rx = 0; rx < Rx; rx++){
				received_vector[rx][0] = DFT_signal[index][rx][f_index][0];
				received_vector[rx][1] = DFT_signal[index][rx][f_index][1];
			}
			DDWE_Update(received_vector, frequency_signal[f_index], W[f_index][0], W[f_index][0], H[f_index], e, f_index);
		}
#endif // DDWE
	}
#endif	    			
} // receiver()


void DFT(int l, int data_number, double(*received_signal)[BURST][2], double(*DFT_signal)[N][2]){
	int	f_index, time, time_index, carrier_index;
	double time_signal[Nf][2], frequency_signal[Nf][2], symbol[N][2];

	for (time_index = 0; time_index < Nf; time_index++) {
		time = data_number * SAMPLEN + GI + time_index;
		time_signal[time_index][0] = received_signal[l][time][0];
		time_signal[time_index][1] = received_signal[l][time][1];
	}

	FFT(time_signal, frequency_signal);

	for (carrier_index = 0; carrier_index < N; carrier_index++) {
		symbol[carrier_index][0] = frequency_signal[(carrier_index + Nf - N / 2) % Nf][0];
		symbol[carrier_index][1] = frequency_signal[(carrier_index + Nf - N / 2) % Nf][1];
	}

	for (f_index = 0; f_index < N; f_index++){
		DFT_signal[l][f_index][0] = symbol[f_index][0];
		DFT_signal[l][f_index][1] = symbol[f_index][1];
	}
}