#include "const.h"

extern double transmitted_sigmat[TS][BS][BURST][2];
extern double wh_sequence[BS][Np][N], channel_power[N][2], her[N][TS][2];
double h_dsrd[N][TS][2];
double H[N][TS][BS][2];
extern errno_t err;
extern double CNR, EbNo;
extern int loop;
double OutPflag;

void receiver(double(*Received_Signal)[BURST][2], int(*rbit)[Nd][N_dbps], double(*wh_sequence)[Np][N]){
	FILE *mse;
	char strI[10];
	char File_Name[20] = "mse_";
	char Ext[5] = ".dat";
	int index, ts, k, f_index;
	int coded_bit[KI][N_cbps];
	double received_vector[TS][2], DFT_signal[PACKETN][TS][N][2], e[2];
	double DFT_signal_ideal[PACKETN][TS][N][2];	/* Ideal one is for RLS. */
	double W[N][BS][CL][2];
	double e_rls[N][NUMBER_OF_STEPS][2], e_rls_av[NUMBER_OF_STEPS];
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
			DFT(ts, index, Received_Signal, DFT_signal[index]);
		}
		DFT(0, index, transmitted_sigmat[0], DFT_signal_ideal[index]);
	}
#endif // 			
	RLS_initial_weigth(0, DFT_signal, DFT_signal_ideal, W, H, wh_sequence, e_rls);



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
			DFT(ts, index, Received_Signal, DFT_signal[index]);
		}
		CoherentDetection4HD(index, DFT_signal[index], W, diLLR, Decision);
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
			}
			DDWE_Update(
				received_vector,
				frequency_signal_hd[0][f_index],
				W[f_index][0],
				W[f_index][0],
				H[f_index],
				e_rls,
				frequency_signal_hd[1][f_index],
				f_index,
				index);
		}
#endif // DDWE
	}				  
	 _itoa(loop, strI, 10);
	strcat(File_Name,strI);
	strcat(File_Name, Ext);
	if ((err = fopen_s(&mse, File_Name, "w+")) != 0)	{
		printf("The file mse.dat file was not opened\n");
		getchar();
		exit(1);
	};
	File_Name[4] = '\0';
	//printf("%d\t%d\n",loop,EbNo);
	for (size_t step = 0; step < Np + Nd; step++)
	{
		e_rls_av[step] = 0;
		for (size_t f_index = 0; f_index < N; f_index++)
		{
			e_rls_av[step] += sqrt(sqr(e_rls[f_index][step][0]) + sqr(e_rls[f_index][step][1]));
			fprintf(mse, "%d", step);
			fprintf(mse, "\t%1.1g", sqrt(sqr(e_rls[f_index][step][0]) + sqr(e_rls[f_index][step][1])));
		}
		e_rls_av[step] = e_rls_av[step] / N;
		//fprintf(mse, "%d\t%e\n", step, e_rls_av[step]);
		fprintf(mse, "\n");
	}
	fclose(mse);		 
	//if (loop == OutPflag + 1) getchar();   
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