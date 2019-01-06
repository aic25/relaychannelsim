#include "const.h"

extern double transmitted_sigmat[TS][BS][BURST][2];
extern double wh_sequence[BS][Np][N][2], channel_power[N][2], her[N][TS][2];
double h_dsrd[N][TS][2];
double H[N][TS][BS][2];
extern errno_t err;
extern double CNR, EbNo;
extern int loop;
double OutPflag;
extern double e_rls[N][NUMBER_OF_STEPS][2], e_rls_avg[NUMBER_OF_STEPS];
extern double e_w[N][NUMBER_OF_STEPS][2], e_w_avg[NUMBER_OF_STEPS];
extern double e_w3[N][NUMBER_OF_STEPS][2], e_w3_avg[NUMBER_OF_STEPS];

void receiver(
	int(*R_index),
	double(*Received_Signal)[BURST][2],
	double(*Received_Signal_wonoise)[BURST][2],
	int(*rbit)[Nd][N_dbps],
	double(*wh_sequence)[Np][N][2]
	){

	FILE *mse, *w, *w3, *corr_eVw, *corr_eVw3;
	char strI[10];
	char File_Name[15] = "mse_";
	char File_Name_forW[15] = "W_";
	char File_Name_forW3[15] = "W3_";
	char File_Name_for_corr_eVW[25] = "corr_eVw_";
	char File_Name_for_corr_eVW3[25] = "corr_eVw3_";
	char File_Name_for_corr_WVW3[25] = "corr_wVw3_";
	char Ext[5] = ".dat";
	int index, ts, k, f_index, Rev_op, Cor_op;
	int coded_bit[KI][N_cbps];
	double received_vector[TS][2], DFT_signal[PACKETN][TS][N][2], e[2];
	double received_vector_wonoise[TS][2], DFT_signal_wonoise[PACKETN][TS][N][2];

	double DFT_signal_ideal[PACKETN][TS][N][2];	/* Ideal one is for RLS. */
	double W[N][BS][CL][2];
	double temp[2], temp1[2];
	int Decision[KI];

#ifdef SD
	double diLLR[Nd][N_cbps];	/* Pre-estimation buffer? */
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
	channel_estimation4ZC_4all_ch(Received_Signal, H);
#else	/* For the case of WH */
	for (index = 0; index < Np; index++){
		for (ts = 0; ts < TS; ts++) {
			DFT(ts, index, Received_Signal, DFT_signal[index]);
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
		for (ts = 0; ts < TS; ts++) {
#ifdef NOISE
			DFT(ts, index, Received_Signal, DFT_signal[index]);
#else
			DFT(ts, index, Received_Signal_wonoise, DFT_signal[index]);
#endif
		}
		DFT(0, index, transmitted_sigmat[1], DFT_signal_ideal[index]);
	}
#endif //
#ifdef W_SELECTIVE 								
	Rev_op = 0; Cor_op = 0;
#ifdef W3
	for (f_index = 0; f_index < N; f_index++)
	{
		RLS_initial_weigth(f_index, R_index[f_index], DFT_signal, DFT_signal_ideal, W, H, wh_sequence, e_rls, e_w, e_w3);

		ScalarMulScalar_H(W[f_index][R_index[f_index]][2], W[f_index][R_index[f_index]][2], temp);
		if (sqr(temp[0]) + sqr(temp[1]) > CORRECTION_TRESHOLD)
		{
			Rev_op++; R_index[f_index] = 1;
			//cout << "First vector power: " << sqr(temp[0]) + sqr(temp[1]) << endl;
			RLS_initial_weigth(f_index, R_index[f_index], DFT_signal, DFT_signal_ideal, W, H, wh_sequence, e_rls, e_w, e_w3);
			ScalarMulScalar_H(W[f_index][R_index[f_index]][2], W[f_index][R_index[f_index]][2], temp1);
			//	cout << "Second vector power: " << sqr(temp1[0]) + sqr(temp1[1]) << endl;
			if ((sqr(temp1[0]) + sqr(temp1[1])) > (sqr(temp[0]) + sqr(temp[1])))
			{
				Cor_op++; R_index[f_index] = 0;
				RLS_initial_weigth(f_index, R_index[f_index], DFT_signal, DFT_signal_ideal, W, H, wh_sequence, e_rls, e_w, e_w3);
			}
			//getchar(); 
		}
	}
#else
	for (f_index = 0; f_index < N; f_index++)
	{
		RLS_initial_weigth(f_index, R_index[f_index], DFT_signal, DFT_signal_ideal, W, H, wh_sequence, e_rls, e_w, e_w3);

		VectorMulVectorToScalarH_CL(W[f_index][R_index[f_index]], W[f_index][R_index[f_index]], temp);
		if (sqr(temp[0]) + sqr(temp[1]) > CORRECTION_TRESHOLD)
		{
			Rev_op++; R_index[f_index] = 1;
			//cout << "First vector power: " << sqr(temp[0]) + sqr(temp[1]) << endl;
			RLS_initial_weigth(f_index, R_index[f_index], DFT_signal, DFT_signal_ideal, W, H, wh_sequence, e_rls, e_w, e_w3);
			VectorMulVectorToScalarH_CL(W[f_index][R_index[f_index]], W[f_index][R_index[f_index]], temp1);
			//	cout << "Second vector power: " << sqr(temp1[0]) + sqr(temp1[1]) << endl;
			if ((sqr(temp1[0]) + sqr(temp1[1])) > (sqr(temp[0]) + sqr(temp[1])))
			{
				Cor_op++; R_index[f_index] = 0;
				RLS_initial_weigth(f_index, R_index[f_index], DFT_signal, DFT_signal_ideal, W, H, wh_sequence, e_rls, e_w, e_w3);
			}
			//getchar(); 
		}
	}
#endif // W3

	cout << "#REV: " << Rev_op << ",\t#CORR: " << Cor_op << endl;	 
#else
	for (f_index = 0; f_index < N; f_index++)
	{
		RLS_initial_weigth(f_index, R_index[f_index], DFT_signal, DFT_signal_ideal, W, H, wh_sequence, e_rls, e_w, e_w3);
	}
#endif // W_SELECTIVE 
#endif	// RLS

	/****************************** DETECTION ******************************/
#ifdef SD
	double frequency_signal_sd[N][2];
	for (index = Np; index < PACKETN; index++){
		for (ts = 0; ts < TS; ts++) {
			DFT(ts, index, Received_Signal, DFT_signal[index]);
		}
		CoherentDetection4SD(index, DFT_signal[index], W, diLLR);
		k = index - Np;
		SD_Viterbi_decoder(diLLR[k], rbit[k]);

#ifdef DDWE
		convolutional_encoder(rbit[k], coded_bit);	/* k is the symbol index */
		interleaver(coded_bit, coded_bit);
		data_insert_regenerated(k, coded_bit, frequency_signal_sd);
		for (f_index = 0; f_index < N; f_index++){
			for (ts = 0; ts < TS; ts++){
				received_vector[ts][0] = DFT_signal[index][ts][f_index][0];
				received_vector[ts][1] = DFT_signal[index][ts][f_index][1];
			}
			DDWE_Update(received_vector, frequency_signal_sd[f_index], W[f_index][0], W[f_index][0], H[f_index], e, f_index);
		}
#endif // DDWE
	}
#endif

#ifdef HD
	double frequency_signal_hd[KI][N][2];
	for (index = Np; index < PACKETN; index++){
		for (ts = 0; ts < TS; ts++) {											
			DFT(ts, index, Received_Signal_wonoise, DFT_signal_wonoise[index]);
#ifdef NOISE
			DFT(ts, index, Received_Signal, DFT_signal[index]);
#else
			DFT(ts, index, Received_Signal_wonoise, DFT_signal[index]);
#endif
		}
		CoherentDetection4HD(R_index, index, DFT_signal[index], W, diLLR, Decision);
		k = index - Np;										 /* k is the symbol index */
		HD_Viterbi_decoder(diLLR[0][k], rbit[0][k]);
		HD_Viterbi_decoder(diLLR[1][k], rbit[1][k]);

#ifdef DDWE												
		convolutional_encoder(rbit[0][k], coded_bit[0]);
		convolutional_encoder(rbit[1][k], coded_bit[1]);
		interleaver(coded_bit[0], coded_bit[0]);
		interleaver(coded_bit[1], coded_bit[1]);
		data_insert_regenerated(k, coded_bit[0], frequency_signal_hd[0]);
		data_insert_regenerated(k, coded_bit[1], frequency_signal_hd[1]);
		for (f_index = 0; f_index < N; f_index++){
			for (ts = 0; ts < TS; ts++){
				received_vector[ts][0] = DFT_signal[index][ts][f_index][0];
				received_vector[ts][1] = DFT_signal[index][ts][f_index][1];
				received_vector_wonoise[ts][0] = DFT_signal_wonoise[index][ts][f_index][0];
				received_vector_wonoise[ts][1] = DFT_signal_wonoise[index][ts][f_index][1];
			}
			DDWE_Update(
				Decision,
				R_index[f_index],
				received_vector,
				received_vector_wonoise,
				frequency_signal_hd[R_index[f_index]][f_index],
				W[f_index][R_index[f_index]],
				W[f_index][R_index[f_index]],
				H[f_index],
				e_rls,
				e_w,
				e_w3,
				frequency_signal_hd[(R_index[f_index] + 1) % 2][f_index],
				f_index,
				index);
		}
#endif // DDWE
	}

#ifdef ERR_POUT_PLOOP
	Spatio_Temporal_Printout(loop, e_rls_avg, e_rls, File_Name);
	Spatio_Temporal_Printout(loop, e_w_avg, e_w, File_Name_forW);
	Spatio_Temporal_Printout(loop, e_w3_avg, e_w3, File_Name_forW3);
#endif
#ifdef CORR_POUT_PLOOP
	Correlation(EbNo, loop, e_rls, e_w, File_Name_for_corr_eVW);
	Correlation(EbNo, loop, e_rls, e_w3, File_Name_for_corr_eVW3);
	Correlation(EbNo, loop, e_w, e_w3, File_Name_for_corr_WVW3);
#endif // CORR_POUT_PLOOP							  

#endif	    			
} // receiver()


void DFT(int l, int data_number, double(*Received_Signal)[BURST][2], double(*DFT_signal)[N][2]){
	int	f_index, time, time_index, carrier_index;
	double time_signal[Nf][2], frequency_signal[Nf][2], symbol[N][2];

	for (time_index = 0; time_index < Nf; time_index++) {
		time = data_number * SAMPLEN + GI + time_index;
		time_signal[time_index][0] = Received_Signal[l][time][0];
		time_signal[time_index][1] = Received_Signal[l][time][1];
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