#include "const.h"
void transmitter(int(*bit)[N_dbps], double(*transmitted_signal)[2], int bs, double(*BS_signal_power)){

	int k, t_index, time, nd;
	int coded_bit[Nd][N_cbps];
	double time_signal[Nd][Nf][2], data_signal[SAMPLEN][2];
	double frequency_signal[Nd][Nf][2];

	BS_signal_power[bs] = 0;

	/* Those for loops look like they can be combined into one */
	for (k = 0; k < Nd; k++){
		info_generator(bit[k]);	/* generate one and zeros */
		convolutional_encoder(bit[k], coded_bit[k]);	/* k is the symbol index */
		interleaver(coded_bit[k], coded_bit[k]);
	}

	for (nd = 0; nd < Nd; nd++){
		data_insert(nd, coded_bit, frequency_signal);	/* Modulate bits to symbols for each packet */
		IFFT(frequency_signal[nd], time_signal[nd]);	/* Generate OFDM symbol of index nd */
	}

	for (nd = 0; nd < Nd; nd++){
		data_GI_insert(time_signal[nd], data_signal);
		for (t_index = 0; t_index < SAMPLEN; t_index++){
			time = (Np + nd) * SAMPLEN + t_index;
			transmitted_signal[time][0] = data_signal[t_index][0];
			transmitted_signal[time][1] = data_signal[t_index][1];
		}
	}
	for ( k = 0; k < BURST; k++){
		BS_signal_power[bs] += sqr(transmitted_signal[k][0]) + sqr(transmitted_signal[k][1]);
	}
	BS_signal_power[bs] /= BURST;
}

void info_generator(int(*data_bit)){
	int n, i, number;

	for (n = 0; n < N_dbps - TAIL;){
		number = rand();
		for (i = 0; i < 15 && n < (N_dbps - TAIL); i++, n++) data_bit[n] = (number >> i) & 0x1;	/* Why generate this way? */
	}
	for (n = 0; n < TAIL; n++) data_bit[N_dbps - TAIL + n] = 0;
}

void data_insert(int nd, int(*bit)[N_cbps], double(*frequency_signal)[Nf][2]){
	int n, index;
	double bin2sgnlQPSK[2] = { -1.0, 1.0 };
	double bin2sgnlQAM[2][2] = { { -3.0, -1.0 }, { 3.0, 1.0 } };
	double bin2sgnl64QAM[2][2][2] = {
		{ { -7.0, -5.0 }, { -1.0, -3.0 } },
		{ { 7.0, 5.0 }, { 1.0, 3.0 } },
	};
	double transmitted_symbol[N][2];
	double FFT_symbol[N][2];

	index = 0;
	for (n = 0; n < N; n++) {
		if (MODULATION == QPSK){//QPSK
			transmitted_symbol[n][0] = OneBySqrt2 * bin2sgnlQPSK[bit[nd][index]];
			transmitted_symbol[n][1] = OneBySqrt2 * bin2sgnlQPSK[bit[nd][index + 1]];
		}
		else if (MODULATION == QAM16){//16QAM
			transmitted_symbol[n][0] = OneBySqrt10 * bin2sgnlQAM[bit[nd][index]][bit[nd][index + 1]];
			transmitted_symbol[n][1] = OneBySqrt10 * bin2sgnlQAM[bit[nd][index + 2]][bit[nd][index + 3]];
		}
		else if (MODULATION == QAM64){//64QAM
			transmitted_symbol[n][0] = OneBySqrt42 * bin2sgnl64QAM[bit[nd][index]][bit[nd][index + 1]][bit[nd][index + 2]];
			transmitted_symbol[n][1] = OneBySqrt42 * bin2sgnl64QAM[bit[nd][index + 3]][bit[nd][index + 4]][bit[nd][index + 5]];
		}
		index += MODULATION;
	}

	for (n = 0; n < N; n++) {
		FFT_symbol[n][0] = transmitted_symbol[n][0];
		FFT_symbol[n][1] = transmitted_symbol[n][1];
	}

	for (n = 0; n < Nf; n++) frequency_signal[nd][n][0] = frequency_signal[nd][n][1] = 0.0;

	for (n = 0; n < N; n++) {
		frequency_signal[nd][(n + Nf - N / 2) % Nf][0] = FFT_symbol[n][0];
		frequency_signal[nd][(n + Nf - N / 2) % Nf][1] = FFT_symbol[n][1];
	}
}

void data_insert_regenerated(int nd, int(*bit), double(*frequency_signal_out)[2]){
	int n, index, carrier_index;
	double bin2sgnlQPSK[2] = { -1.0, 1.0 };
	double bin2sgnlQAM[2][2] = { { -3.0, -1.0 }, { 3.0, 1.0 } };
	double bin2sgnl64QAM[2][2][2] = {
		{ { -7.0, -5.0 }, { -1.0, -3.0 } },
		{ { 7.0, 5.0 }, { 1.0, 3.0 } },
	};
	double transmitted_symbol[N][2], frequency_signal[Nf][2];
	double FFT_symbol[N][2];

	index = 0;
	for (n = 0; n < N; n++) {
		if (MODULATION == QPSK){//QPSK
			transmitted_symbol[n][0] = OneBySqrt2 * bin2sgnlQPSK[bit[index]];
			transmitted_symbol[n][1] = OneBySqrt2 * bin2sgnlQPSK[bit[index + 1]];
		}
		else if (MODULATION == QAM16){//16QAM
			transmitted_symbol[n][0] = OneBySqrt10 * bin2sgnlQAM[bit[index]][bit[index + 1]];
			transmitted_symbol[n][1] = OneBySqrt10 * bin2sgnlQAM[bit[index + 2]][bit[index + 3]];
		}
		else if (MODULATION == QAM64){//64QAM
			transmitted_symbol[n][0] = OneBySqrt42 * bin2sgnl64QAM[bit[index]][bit[index + 1]][bit[index + 2]];
			transmitted_symbol[n][1] = OneBySqrt42 * bin2sgnl64QAM[bit[index + 3]][bit[index + 4]][bit[index + 5]];
		}
		index += MODULATION;
	}

	for (n = 0; n < N; n++) {
		FFT_symbol[n][0] = transmitted_symbol[n][0];
		FFT_symbol[n][1] = transmitted_symbol[n][1];
	}

	for (n = 0; n < Nf; n++) frequency_signal[n][0] = frequency_signal[n][1] = 0.0;

	for (n = 0; n < N; n++) {
		frequency_signal[(n + Nf - N / 2) % Nf][0] = FFT_symbol[n][0];
		frequency_signal[(n + Nf - N / 2) % Nf][1] = FFT_symbol[n][1];
	}
	for ( carrier_index = 0; carrier_index < N; carrier_index++) {
		frequency_signal_out[carrier_index][0] = frequency_signal[(carrier_index + Nf - N / 2) % Nf][0];
		frequency_signal_out[carrier_index][1] = frequency_signal[(carrier_index + Nf - N / 2) % Nf][1];
	}
}

void data_GI_insert(double(*time_signal)[2], double(*out_signal)[2]){
	int k;

	for (k = 0; k < GI; k++){
		out_signal[k][0] = time_signal[Nf - GI + k][0];
		out_signal[k][1] = time_signal[Nf - GI + k][1];	
		/*out_signal[k][0] = 0;
		out_signal[k][1] = 0;	 */
	}
	for (k = GI; k < SAMPLEN; k++){
		out_signal[k][0] = time_signal[k - GI][0];
		out_signal[k][1] = time_signal[k - GI][1];
	}
}


