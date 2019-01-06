#include	"const.h"
void Relay_transmitter(int(*bit)[N_dbps], double(*transmitted_signal)[2]){
	int k, t_index, time, nd;
	int coded_bit[Nd][N_cbps];
	double frequency_signal[Nd][Nf][2], time_signal[Nd][Nf][2], data_signal[SAMPLEN][2];

	for (k = 0; k < Nd; k++){
		convolutional_encoder(bit[k], coded_bit[k]);
		interleaver(coded_bit[k], coded_bit[k]);
	}

	for (nd = 0; nd < Nd; nd++){
		data_insert(nd, coded_bit, frequency_signal);
		IFFT(frequency_signal[nd], time_signal[nd]);
	}

	for (nd = 0; nd < Nd; nd++){
		data_GI_insert(time_signal[nd], data_signal);
		for (t_index = 0; t_index < SAMPLEN; t_index++){
			time = (Np + nd) * SAMPLEN + t_index;
			transmitted_signal[time][0] = data_signal[t_index][0];
			transmitted_signal[time][1] = data_signal[t_index][1];
		}
	}
}
