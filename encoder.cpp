#include	"const-function-declerations.h"		

int HSIZE, VSIZE;	/* Horizontal SIZE and Vertical SIZE */

void convolutional_encoder(int(*bit), int(*code)){
	int n, k, s[6], index, state = 0;
	int temp_code[2 * N_dbps];

	for (index = 0; index < N_dbps; index++){
		for (k = 0; k < 6; k++) s[k] = (state >> k) & 0x1;
		/* encoder1 */
		temp_code[2 * index] = bit[index] ^ s[4] ^ s[3] ^ s[1] ^ s[0];
		/* encoder2 */
		temp_code[2 * index + 1] = bit[index] ^ s[5] ^ s[4] ^ s[3] ^ s[0];

		state = (bit[index] << 5) + (s[5] << 4) + (s[4] << 3) + (s[3] << 2) + (s[2] << 1) + s[1];
	}
	/********************** no puncturing ********************************/
	for (n = 0; n < 2 * N_dbps; n++){
		code[n] = temp_code[n];
	}
}

void interleaver(int(*input), int(*output)){
	int store[N_cbps];
	int n, i, temp;
	int size_constant;
	size_constant = (N == 512) ? 32 : ((N == 128) ? 16 : 8);

	if (MODULATION == QPSK){	 // 32-48-64
		HSIZE = size_constant;
		VSIZE = size_constant;
	}
	else if (MODULATION == QAM16){
		HSIZE = size_constant;
		VSIZE = 2 * size_constant;
	}
	else{
		HSIZE = (int)(1.5*size_constant);
		VSIZE = 2 * size_constant;
	}

	for (n = 0; n < N_cbps; n++) store[n] = input[n];
	for (n = 0; n < N_cbps; n++){
		temp = n % HSIZE;
		i = VSIZE*temp + n / HSIZE;
		output[+i] = store[n];
	}
}
