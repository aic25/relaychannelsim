#include "const-function-declerations.h"	 	
extern double transmitted_sigmat[TS][BS][BURST][2];
extern double wh_sequence[BS][Np][N][2], channel_power[N][2], her[N][TS][2];
double h_dsrd[N][TS][2];
double H[N][TS][BS][2];
extern errno_t err;
extern double CNR, EbNo;
extern int loop, R_index[N];

void receiver(double(*received_signal)[BURST][2], int(*rbit)[Nd][N_dbps]){
	int index, ts, j, k, f_index, coded_bit[KI][N_cbps], Decision[KI], Rev_op, Cor_op;
	double received_vector[TS][2], DFT_signal[PACKETN][TS][N][2], temp[2], temp1[2];
	double DFT_signal_ideal[PACKETN][TS][N][2], e_rls[N][NUMBER_OF_STEPS][2];	/* Ideal one is for RLS. */
	double W[N][BS][CL][2];
	double temp_mat[N][TS][BS][2];
	int Swap_UD[N];
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
	channel_estimation4ZC_4all_ch(received_signal, H);
#else	/* For the case of WH */
	for (index = 0; index < Np; index++){
		for (ts = 0; ts < TS; ts++) {
			DFT(ts, index, received_signal, DFT_signal[index]);
		}
		DFT(0, index, transmitted_sigmat[0], DFT_signal_ideal[index]);
	}
	WH_channel_estimation(DFT_signal, H, wh_sequence);
#ifdef LSCHup
	for (f_index = 0; f_index < N; f_index++){
		MatrixEquate(temp_mat[f_index], H[f_index]);
	}
#endif // LSCHup

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
			DFT(ts, index, received_signal, DFT_signal[index]);
		}
		DFT(0, index, transmitted_sigmat[0], DFT_signal_ideal[index]);
	}
#endif // 			
	//RLS_CombCoeff(0, DFT_signal, DFT_signal_ideal, W, H);
	/**/
#ifndef UD_SELECT
	for (f_index = 0; f_index < N; f_index++)
	{
		R_index[f_index] = 0;
		Swap_UD[f_index] = 0;
		RLS_CombCoeff(f_index, R_index[f_index], Swap_UD[f_index], DFT_signal, DFT_signal_ideal, W, H, e_rls);
	}						
#else
	//Rev_op = 0; Cor_op = 0;
	for (f_index = 0; f_index < N; f_index++){
		R_index[f_index] = 0;
		Swap_UD[f_index] = 0;
		RLS_CombCoeff(f_index, R_index[f_index], Swap_UD[f_index], DFT_signal, DFT_signal_ideal, W, H, e_rls, 1);
		RLS_CombCoeff(f_index, R_index[f_index], Swap_UD[f_index], DFT_signal, DFT_signal_ideal, W, H, e_rls, 0);
#ifndef E_SELECT
		ScalarMulScalar_H(W[f_index][R_index[f_index]][0], W[f_index][R_index[f_index]][0], temp);
#else
		//temp[0] = e_rls[f_index][PREAMBLE_LENGTH - 1][0]; temp[1] = e_rls[f_index][PREAMBLE_LENGTH - 1][1];
		VectorMulVectorToScalarH_Np(e_rls[f_index], e_rls[f_index], temp);
#endif
		if (sqr(temp[0]) + sqr(temp[1]) > 0){
			Swap_UD[f_index] = 1;
			RLS_CombCoeff(f_index, R_index[f_index], Swap_UD[f_index], DFT_signal, DFT_signal_ideal, W, H, e_rls, 1);
			RLS_CombCoeff(f_index, R_index[f_index], Swap_UD[f_index], DFT_signal, DFT_signal_ideal, W, H, e_rls, 0);
#ifndef E_SELECT
			ScalarMulScalar_H(W[f_index][R_index[f_index]][0], W[f_index][R_index[f_index]][0], temp1);
#else
			//temp1[0] = e_rls[f_index][PREAMBLE_LENGTH - 1][0]; temp1[1] = e_rls[f_index][PREAMBLE_LENGTH - 1][1];
			VectorMulVectorToScalarH_Np(e_rls[f_index], e_rls[f_index], temp1);
#endif
			if ((sqr(temp1[0]) + sqr(temp1[1])) > (sqr(temp[0]) + sqr(temp[1]))){
				Swap_UD[f_index] = 0;
				RLS_CombCoeff(f_index, R_index[f_index], Swap_UD[f_index], DFT_signal, DFT_signal_ideal, W, H, e_rls, 1);
			}
		}
	}
#endif
#endif	// RLS
	/****************************** DETECTION ******************************/
#ifdef SD
	double frequency_signal[N][2];
	for (index = Np; index < PACKETN; index++){
		for (ts = 0; ts < TS; ts++) {
			DFT(ts, index, received_signal, DFT_signal[index]);
		}
		CoherentDetection4SD(index, DFT_signal[index], W, diLLR);
		k = index - Np;
		SD_Viterbi_decoder(diLLR[k], rbit[k]);

#ifdef DDWE
		convolutional_encoder(rbit[k], coded_bit);	/* k is the symbol index */
		interleaver(coded_bit, coded_bit);
		data_insert_regenerated(k, coded_bit, frequency_signal);
		for (f_index = 0; f_index < N; f_index++){
			for (ts = 0; ts < TS; ts++){
				received_vector[ts][0] = DFT_signal[index][ts][f_index][0];
				received_vector[ts][1] = DFT_signal[index][ts][f_index][1];
			}
			DF_CombCoeff(received_vector, frequency_signal[f_index], W[f_index][0], W[f_index][0], H[f_index], e, f_index);
		}
#endif // DDWE
	}			   
#endif

#ifdef HD
	double frequency_signal[KI][N][2];
	for (index = Np; index < PACKETN; index++){
		for (ts = 0; ts < TS; ts++) {
			DFT(ts, index, received_signal, DFT_signal[index]);
		}
		CoherentDetection4HD(index, DFT_signal[index], W, H, diLLR, Decision, R_index, Swap_UD);
		k = index - Np;
		if (CHANNEL_ENCODING) {
			HD_Viterbi_decoder(diLLR[0][k], rbit[0][k]);
			HD_Viterbi_decoder(diLLR[1][k], rbit[1][k]);
		}
		else {
			for (j = 0; j < N_dbps; j++) {
				rbit[0][k][j] = diLLR[0][k][j];
				rbit[1][k][j] = diLLR[1][k][j];
			}
		}

#ifdef DDWE
		convolutional_encoder(rbit[0][k], coded_bit[0]);
		interleaver(coded_bit[0], coded_bit[0]);
		data_insert_regenerated(k, coded_bit[0], frequency_signal[0]);

		convolutional_encoder(rbit[1][k], coded_bit[1]);
		interleaver(coded_bit[1], coded_bit[1]);
		data_insert_regenerated(k, coded_bit[1], frequency_signal[1]);

		for (f_index = 0; f_index < N; f_index++){
			for (ts = 0; ts < TS; ts++){
				received_vector[ts][0] = DFT_signal[index][ts][f_index][0];
				received_vector[ts][1] = DFT_signal[index][ts][f_index][1];
			}
#ifdef LSCHup
			//MatrixEquate(temp_mat, H[f_index]);
			WH_chest_update(received_vector, temp_mat[f_index], frequency_signal[R_index[f_index]][f_index], f_index);
			//if (index > Np + 10){ MatrixEquate(H[f_index], temp_mat[f_index]); }
#endif // Wh  
			DF_CombCoeff(index, f_index, Decision, R_index[f_index], Swap_UD[f_index], received_vector,
				frequency_signal[R_index[f_index]][f_index],
				frequency_signal[(R_index[f_index] + 1) % 2][f_index],
				W[f_index][R_index[f_index]], W[f_index][R_index[f_index]], H[f_index], e_rls[f_index][index]);

		}
#endif // DDWE
	}
	RLS_Error(e_rls);
#endif	    			
} // receiver()


void DFT(int l, int data_number, double(*received_signal)[BURST][2], double(*DFT_signal)[N][2]){
	int	f_index, time, time_index, carrier_index;
	double time_signal[Nf][2], frequency_signal[Nf][2], symbol[N][2];

	for (time_index = 0; time_index < Nf; time_index++){
		time = data_number * SAMPLEN + GI + time_index;
		time_signal[time_index][0] = received_signal[l][time][0];
		time_signal[time_index][1] = received_signal[l][time][1];
	}

	(Fourier_MOD == ON) ? FFT(time_signal, frequency_signal) : DFT(time_signal, frequency_signal);

	for (carrier_index = 0; carrier_index < N; carrier_index++){
		symbol[carrier_index][0] = frequency_signal[(carrier_index + Nf - N / 2) % Nf][0];
		symbol[carrier_index][1] = frequency_signal[(carrier_index + Nf - N / 2) % Nf][1];
	}

	for (f_index = 0; f_index < N; f_index++){
		DFT_signal[l][f_index][0] = symbol[f_index][0];
		DFT_signal[l][f_index][1] = symbol[f_index][1];
	}
}

void RLS_Error(double e_rls[N][NUMBER_OF_STEPS][2]){
	ofstream myfile("MSE/" + std::to_string(loop) + "_mse.dat");
	double e[NUMBER_OF_STEPS];
	int i, j;
	for (i = 0; i < NUMBER_OF_STEPS; i++){
		e[i] = 0;
		for (j = 0; j < N; j++){
			e[i] += (sqr(e_rls[j][i][0]) + sqr(e_rls[j][i][1])) / N;
		}
		myfile << i << "\t" << e[i] << endl;
	}
	myfile.close();
}