#include "const.h"
extern double transmitted_sigmat[TS][BS][BURST][2];
extern double wh_sequence[BS][Np][N], channel_power[N][2], hher[N][TS][2];
double h_dsrd[N][TS][2];
double H[N][CL][BS][2];

void receiver(double(*received_signal)[BURST][2], int(*rbit)[Nd][N_dbps]){
	int np, ts, k, f_index;
	int coded_bit[BS][N_cbps];
	double received_vector[CL][2], DFT_signal[PACKETN][TS][N][2], e[2], 
		s_qpsk[4][2] = { { 1.0, 1.0 }, { 1.0, -1.0 }, { -1.0, -1.0 }, { -1.0, 1.0 } }/*qpsk*/;
	double DFT_signal_ideal[PACKETN][KI][N][2], DFT_signal_ideal_bs3[PACKETN][KI][N][2];	/* Ideal one is for RLS. */
	double W[N][KI][CL][2], tempmat[CL][CL][2];
#ifdef SD
	double diLLR[KI][Nd][N_cbps];	/* Pre-estimation buffer? */
#endif
#ifdef HD
	int diLLR[KI][Nd][N_cbps];	/* Pre-estimation buffer? */
#endif 			  

	/************************** CHANNEL ESTIMATION *************************/
#ifdef CUNNING
	ideal_channel_estimation4all_ch(H);
#endif
#ifndef CUNNING
#ifdef ZC
	channel_estimation4ZC_4all_ch(received_signal, H);
#else	/* For the case of WH */
	for (np = 0; np < Np; np++){
		for (ts = 0; ts < TS; ts++) {
			DFT(ts, np, received_signal, DFT_signal[np]);
		}
		DFT(0, np, transmitted_sigmat[0], DFT_signal_ideal[np]);
	}
	WH_channel_estimation(DFT_signal, H, wh_sequence);
#endif  
#endif // !CUNNING
	h_desired(H, h_dsrd);
	VectorHermite(h_dsrd, hher);
	VectorMulVectorToscalar(hher, h_dsrd, channel_power);

	/*********************** Weight Generation *****************************/
#ifdef MMSE
	MMSE_Weight_Matrix(H, W);
#endif
#ifdef MRC
	mrc_weight(H);
#endif
#ifdef RLS	/* more to add */
#ifndef Wh
	for (np = 0; np < Np; np++){
		for (ts = 0; ts < TS; ts++) {
			DFT(ts, np, received_signal, DFT_signal[np]);
		}
		DFT(0, np, transmitted_sigmat[0], DFT_signal_ideal[np]);
		DFT_forusr3(0, np, transmitted_sigmat[0], DFT_signal_ideal_bs3[np]);
		for (size_t f_index = 0; f_index < N; f_index++){
			DFT_signal_ideal[np][1][f_index][0] = DFT_signal_ideal_bs3[np][0][f_index][0];
		}
	}
#endif // 			
	RLS_initial_weigth(0, DFT_signal, DFT_signal_ideal, W, H);	
	//RLS_initial_weigth(1, DFT_signal, DFT_signal_ideal, W, H);
#endif	// RLS
	/****************************** DETECTION ******************************/
#ifdef SD
	double frequency_signal[KI][N][2];
	for (size_t index = Np; index < PACKETN; index++){
		for (ts = 0; ts < TS; ts++) {
			DFT(ts, index, received_signal, DFT_signal[index]);
		}
		CoherentDetection4SD(0, index, DFT_signal[index], W, diLLR);
		k = index - Np;
		SD_Viterbi_decoder(diLLR[0][k], rbit[0][k]);
		SD_Viterbi_decoder(diLLR[1][k], rbit[1][k]); 
#ifdef DDWE
		convolutional_encoder(rbit[0][k], coded_bit[0]);	/* k is the symbol index */
		convolutional_encoder(rbit[1][k], coded_bit[1]);
		interleaver(coded_bit[0], coded_bit[0]);
		interleaver(coded_bit[1], coded_bit[1]);
		data_insert_regenerated(k, coded_bit[0], frequency_signal[0]);
		data_insert_regenerated(k, coded_bit[1], frequency_signal[1]);
		for (f_index = 0; f_index < N; f_index++){
			for (ts = 0; ts < TS; ts++){
				received_vector[ts][0] = DFT_signal[index][ts][f_index][0];
				received_vector[ts][1] = DFT_signal[index][ts][f_index][1];
			}
			received_vector[CL-1][0] = frequency_signal[1][f_index][0];
			received_vector[CL-1][1] = frequency_signal[1][f_index][1];
			//frequency_signal[f_index][0] = s_qpsk[s1_index][0];
			//frequency_signal[f_index][1] = s_qpsk[s1_index][1];
			DDWE_Update(0, received_vector, frequency_signal[0][f_index], W[f_index][0], W[f_index][0], H[f_index], e, f_index);
		}
#endif // DDWE
	}			   
#endif

#ifdef HD
	double frequency_signal[N][2];
	for (size_t index = Np; index < PACKETN; index++){
		for (ts = 0; ts < TS; ts++) {
			DFT(ts, index, received_signal, DFT_signal[index]);
		}
		CoherentDetection4HD(index, DFT_signal[index], W, diLLR[0]);
		k = index - Np;
		HD_Viterbi_decoder(diLLR[0][k], rbit[k]);

#ifdef DDWE
		convolutional_encoder(rbit[k], coded_bit);	/* k is the symbol index */
		interleaver(coded_bit, coded_bit);
		data_insert_regenerated(k, coded_bit, frequency_signal);
		for (f_index = 0; f_index < N; f_index++){
			for (ts = 0; ts < TS; ts++){
				received_vector[ts][0] = DFT_signal[index][ts][f_index][0];
				received_vector[ts][1] = DFT_signal[index][ts][f_index][1];
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

void DFT_forusr3(int l, int data_number, double(*received_signal)[BURST][2], double(*DFT_signal)[N][2]){
	int	f_index, time, time_index, carrier_index;
	double time_signal[Nf][2], frequency_signal[Nf][2], symbol[N][2];

	for (time_index = 0; time_index < Nf; time_index++) {
		time = data_number * SAMPLEN + GI + time_index;
		time_signal[time_index][0] = received_signal[CL-1][time][0];
		time_signal[time_index][1] = received_signal[CL-1][time][1];
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