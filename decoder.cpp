#include "const.h"
struct StateInfo state_table[STATEN][2];
extern int HSIZE, VSIZE;

void state_model_convolutional(){
	int state, input, k, s[6];

	for (state = 0; state < STATEN; state++){
		for (k = 0; k < 6; k++) s[k] = (state >> k) & 0x1;
		for (input = 0; input < 2; input++){
			/* encoder1 */
			state_table[state][input].output0 = input^s[4] ^ s[3] ^ s[1] ^ s[0];
			/* encoder2 */
			state_table[state][input].output1 = input^s[5] ^ s[4] ^ s[3] ^ s[0];
			state_table[state][input].next_state = (input << 5) + (s[5] << 4) + (s[4] << 3) + (s[3] << 2) + (s[2] << 1) + s[1];
			state_table[state][input].pre_state = (s[4] << 5) + (s[3] << 4) + (s[2] << 3) + (s[1] << 2) + (s[0] << 1) + input;
		}
	}
}

void HD_Viterbi_decoder(int(*decode), int(*received_bit)){
	int n, state, path, s[2], input;
	int Pmetric[STATEN], Phistory[N_dbps][STATEN];
	int Bmetric[2], min_path, Pt[STATEN];

	for (state = 0; state < STATEN; state++){
		if (state == 0) Pmetric[state] = 0;
		else Pmetric[state] = 100;
	}

	for (n = 0; n < N_dbps; n++){
		for (state = 0; state < STATEN; state++){
			for (path = 0; path<2; path++){
				s[path] = state_table[state][path].pre_state;
				input = (state >> 5) & 0x1;
				Bmetric[path] = state_table[s[path]][input].output0^decode[2 * n];
				Bmetric[path] += state_table[s[path]][input].output1^decode[2 * n + 1];
				Bmetric[path] += Pmetric[s[path]];
			}
			min_path = 0;
			if (Bmetric[0] > Bmetric[1]) min_path = 1;
			Phistory[n][state] = s[min_path];
			Pt[state] = Bmetric[min_path];
		}

		for (state = 0; state < STATEN; state++){
			Pmetric[state] = Pt[state];
		}
	}

	state = 0;
	for (n = N_dbps - 1; n >= 0; n--){
		received_bit[n] = (state >> 5) & 0x1;
		state = Phistory[n][state];
	}
}


void SD_Viterbi_decoder(double(*decode_signal), int(*received_bit)){
	int n, state, path, s[2], input;
	double signal[2], Pmetric[STATEN], Pt[STATEN], Bmetric[2];
	int Phistory[N_dbps][STATEN], max_path;
	double bin2sgnl[2] = { -1.0, 1.0 };

	for (state = 0; state < STATEN; state++){
		if (state == 0) Pmetric[state] = 0.0;
		else Pmetric[state] = -1.0e6;
	}

	for (n = 0; n < N_dbps; n++){
		for (state = 0; state < STATEN; state++){
			for (path = 0; path < 2; path++){
				s[path] = state_table[state][path].pre_state;
				input = (state >> 5) & 0x1;
				signal[0] = bin2sgnl[state_table[s[path]][input].output0];
				signal[1] = bin2sgnl[state_table[s[path]][input].output1];

				Bmetric[path] = decode_signal[2 * n] * signal[0];
				Bmetric[path] += decode_signal[2 * n + 1] * signal[1];
				Bmetric[path] += Pmetric[s[path]];
			}
			max_path = 0;
			if (Bmetric[0] < Bmetric[1]) max_path = 1;
			Phistory[n][state] = s[max_path];
			Pt[state] = Bmetric[max_path];
		}

		for (state = 0; state < STATEN; state++){
			Pmetric[state] = Pt[state];
		}
	}

	state = 0;
	for (n = N_dbps - 1; n >= 0; n--){
		received_bit[n] = (state >> 5) & 0x1;
		state = Phistory[n][state];
	}
}

void deinterleaver4SD(double(*input), double(*output)){
	int n, i, temp;
	double store[N_cbps];

	for (n = 0; n < N_cbps; n++){
		store[n] = input[n];
	}
	for (n = 0; n < N_cbps; n++){
		temp = n % VSIZE;
		i = HSIZE*temp + n / VSIZE;
		output[i] = store[n];
	}
}

void deinterleaver4HD(int(*input), int(*output)){
	int n, i, temp, store[N_cbps];

	for (n = 0; n < N_cbps; n++){
		store[n] = input[n];
	}
	for (n = 0; n < N_cbps; n++){
		temp = n % VSIZE;
		i = HSIZE*temp + n / VSIZE;
		output[i] = store[n];
	}
}

