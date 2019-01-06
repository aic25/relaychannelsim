#include <stdio.h>
#include <stdlib.h>
#include "const.h"


errno_t err;

void initialization(double(*noise_power), double(*signal_power), double(*interference_power)){
	*noise_power = 0.0;
	*signal_power = 0.0;
	*interference_power = 0.0;
}

void transmitters_initialize(double(*transmitted_signal)[BURST][2]){
	int i, j;
	for ( j = 0; j < BS; j++){
		for (i = 0; i < BURST; i++){
			transmitted_signal[j][i][0] = transmitted_signal[j][i][1] = 0.0;
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////
void transmitter_initialize(double(*transmitted_signal)[BS][BURST][2]){
	int bs, ts, i;

	for (ts = 0; ts < TS; ts++){
		for (bs = 0; bs < BS; bs++) {
			for (i = 0; i < BURST; i++) {
				transmitted_signal[ts][bs][i][0] = transmitted_signal[ts][bs][i][1] = 0.0;
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////


void setting_print(){
	FILE *fdata;
	if ((err = fopen_s(&fdata, BERFILE, "a")) != 0)	{
#ifdef DBG
		printf("The file ");
		printf(BERFILE);
		printf(" was not opened\n");
#endif // DBG		   
		getchar();
		exit(1);
	};	/* end :: fdata */

	fprintf(fdata, "# path: %1d (equal power), delay_dif: %2d, max_delay: %2d\n", PATH, DELAY, (PATH - 1)*DELAY);
	printf("# path: %1d (equal power), delay_dif: %2d, max_delay: %2d\n", PATH, DELAY, (PATH - 1)*DELAY);

	if (MODULATION == QPSK){
		fprintf(fdata, "# mod: QPSK, R=%2.1f, K=7, BS: %1d, TS: %1d\n", CODING_RATE, BS, TS);
		printf("# mod: QPSK, R=%2.1f, K=7, BS: %1d, TS: %1d\n", CODING_RATE, BS, TS);
	}
	else if (MODULATION == QAM16){
		fprintf(fdata, "# mod: QAM16, R=%2.1f, K=7, BS: %1d, TS: %1d\n", CODING_RATE, BS, TS);
		printf("# mod: QAM16, R=%2.1f, K=7, BS: %1d, TS: %1d\n", CODING_RATE, BS, TS);
	}
	else if (MODULATION == QAM64){
		fprintf(fdata, "# mod: QAM64, R=%2.1f, K=7, BS: %1d, TS: %1d\n", CODING_RATE, BS, TS);
		printf("# mod: QAM64, R=%2.1f, K=7, BS: %1d, TS: %1d\n", CODING_RATE, BS, TS);
	}


#ifdef CUNNING
	fprintf(fdata,"# channel estimation: cunning\n");
	printf("# channel estimation: cunning\n");
#else
#ifdef WH_CHEST 	
	fprintf(fdata, "# channel estimation: WH sequence\n");
	printf("# channel estimation: WH sequence\n");
#else 
	fprintf(fdata, "# channel estimation: ZC sequence\n");
	printf("# channel estimation: ZC sequence\n");
#endif
#endif

#ifndef WH_PREAMBLE
	fprintf(fdata, "# Preamble: ZC sequence\n");
	printf("# Preamble: ZC sequence\n");
#else	
	fprintf(fdata, "# Preamble: WH sequence\n");
	printf("# Preamble: WH sequence\n");
#endif

	//fprintf(fdata, "# Relay SNR: Variable\n");
	//printf("# Relay SNR: Variable\n");

	fprintf(fdata, "# Relay SNR: Fixed to %d dB \n", RLY_CNR);
	printf("# Relay SNR: Fixed to %d dB \n", RLY_CNR);
										 
#ifdef DF
	fprintf(fdata, "# Relay Signalling: DF\n");
	printf("# Relay Signalling: DF\n");
#endif
#ifdef AF
	fprintf(fdata, "# Relay Signalling: AF\n");
	printf("# Relay Signalling: AF\n");
#endif					
#if(INTERFERENCE==ON)  
	fprintf(fdata, "# Interference: ON\n");
	printf("# Interference: ON\n");
#else 
	fprintf(fdata,"# Interference: OFF\n");
	printf("# Interference: OFF\n");
#endif 	 
#if(SWITCH==ON)
	fprintf(fdata, "# Switch: ON\n");
	printf("# Switch: ON\n");
#else
	fprintf(fdata,"# Switch: OFF\n");
	printf("# Switch: OFF\n");
#endif
#ifdef NOISE 
	fprintf(fdata, "# Noise: ON\n");
	printf("# Noise: ON\n");
#else 
	fprintf(fdata,"# Noise: OFF\n");
	printf("# Noise: OFF\n");
#endif 		   
#ifdef MRC 
#ifdef HD
	fprintf(fdata, "# Combining: MRC\n");
	printf("# Combining: MRC\n");
#else 
	fprintf(fdata, "# Error: Use HD for MRC\n");
	printf("# Error: Use HD for MRC\n");
	_getch();
#endif
#endif	   
#ifdef MMSE
	fprintf(fdata, "# Combining: MMSE\n");
	printf("# Combining: MMSE\n");
#endif		
#ifdef RLS
#ifndef DDWE
	fprintf(fdata, "# Combining: RLS\n");
	printf("# Combining: RLS\n");
#else
	fprintf(fdata, "# Combining: DDWE\n");
	printf("# Combining: DDWE\n");
#endif // !DDWE	
#endif	 
#ifdef SD
	fprintf(fdata, "# CD: SDVD\n");
	printf("# CD: SDVD\n");
#endif
#ifdef HD
	fprintf(fdata, "# CD: HDVD\n");
	printf("# CD: HDVD\n");
#endif		  
	fprintf(fdata, "# Power loss $\\ksi$: %d dB \n", LOSS);
	printf("# Power loss $\\ksi$: %d dB \n", LOSS);
	fprintf(fdata, "# Preamble length: %d\n", PREAMBLE_LENGTH);
	printf("# Preamble length: %d\n", PREAMBLE_LENGTH);
	fprintf(fdata, "# Initial Eb/No: %d\n", EbNo_FROM);
	printf("# Initial Eb/No: %d\n", EbNo_FROM);
#ifdef W_SELECTIVE

	fprintf(fdata, "# Error function selection: Active\n");
	printf("# Error function selection: Active\n");
#ifdef W3
	fprintf(fdata, "# Selection basis: W3\n");
	printf("# Selection basis: W3\n");
#else
	fprintf(fdata, "# Selection basis: W\n");
	printf("# Selection basis: W\n");
#endif
	fprintf(fdata, "# Correction treshold: %f\n", CORRECTION_TRESHOLD);
	printf("# Correction treshold: %f\n\n", CORRECTION_TRESHOLD);
#endif // W_SELECTIVE



	fclose(fdata);
}
