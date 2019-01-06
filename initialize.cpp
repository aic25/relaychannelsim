#include "const-function-declerations.h"	  

double noise_power, signal_power;
errno_t err;
double sum_rs, path_rs[RS_PATH][RS_RAY_COMPONENT], sum_ms, path_ms[PATH][RAY_COMPONENT];

void initialization(void){
	noise_power = 0.0;
	signal_power = 0.0;
}

void init_channel_models(){

	int path_index, component_index;
	/**/
	sum_rs = 0.0;
	for (path_index = 0; path_index < RS_PATH; path_index++){
		for (component_index = 0; component_index < (path_index == 0 ? (RLY_LOS == ON ? 1 : RS_RAY_COMPONENT) : RS_RAY_COMPONENT); component_index++){
			path_rs[path_index][component_index] = pow(CONST_E, (-path_index / (double)RS_TAU));
			//path_rs[path_index][component_index] = ((RLY_LOS == ON) && (component_index == 0) && (path_index == 0)) ? path_rs[path_index][component_index] : (path_rs[path_index][component_index] / (pow(10, COMPONENT_K / 10)*((double)RS_RAY_COMPONENT - 1)));
			path_rs[path_index][component_index] = ((RLY_LOS == ON) && (component_index == 0) && (path_index == 0)) ? path_rs[path_index][component_index] : (path_rs[path_index][component_index] / RS_RAY_COMPONENT);
			sum_rs += path_rs[path_index][component_index];
		}
	}
	for (path_index = 0; path_index < RS_PATH; path_index++){
		for (component_index = 0; component_index < (path_index == 0 ? (RLY_LOS == ON ? 1 : RS_RAY_COMPONENT) : RS_RAY_COMPONENT); component_index++){
			path_rs[path_index][component_index] /= sum_rs;
		}
	}

	sum_ms = 0.0;
	for (path_index = 0; path_index < PATH; path_index++){
		for (component_index = 0; component_index < (path_index == 0 ? (LOS == ON ? 1 : RAY_COMPONENT) : RAY_COMPONENT); component_index++){
			path_ms[path_index][component_index] = pow(CONST_E, (-path_index / (double)TAU));
			path_ms[path_index][component_index] = ((LOS == ON) && (component_index == 0) && (path_index == 0)) ? path_ms[path_index][component_index] : (path_ms[path_index][component_index] / RAY_COMPONENT);
			sum_ms += path_ms[path_index][component_index];
		}
	}
	for (path_index = 0; path_index < PATH; path_index++){
		for (component_index = 0; component_index < (path_index == 0 ? (LOS == ON ? 1 : RAY_COMPONENT) : RAY_COMPONENT); component_index++){
			path_ms[path_index][component_index] /= sum_ms;
		}
	}
}

void transmitters_initialize(double(*transmitted_signal)[BURST][2]){
	int i;
	for (size_t j = 0; j < BS; j++){
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
		printf("The file ");
		printf(BERFILE);
		printf(" was not opened\n");
		getchar();
		exit(1);
	};	/* end :: fdata */

	fprintf(fdata, "# MS path: %1d (equal power), delay_dif: %2d, max_delay: %2d\n", PATH, DELAY, (PATH - 1)*DELAY);
	printf("# MS path: %1d (equal power), delay_dif: %2d, max_delay: %2d\n", PATH, DELAY, (PATH - 1)*DELAY);
	fprintf(fdata, "# RS path: %1d (equal power), delay_dif: %2d, max_delay: %2d\n", RS_PATH, RS_DELAY, (RS_PATH - 1)*RS_DELAY);
	printf("# RS path: %1d (equal power), delay_dif: %2d, max_delay: %2d\n", RS_PATH, RS_DELAY, (RS_PATH - 1)*RS_DELAY);
	fprintf(fdata, "Np=%3d, PL= %3d, LOSS=%3d, RLY_CNR=%2d\n", Np, PREAMBLE_LENGTH, LOSS, RLY_CNR);
	printf("Np=%3d, PL= %3d, LOSS=%3d, RLY_CNR=%2d\n", Np, PREAMBLE_LENGTH, LOSS, RLY_CNR);

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
	fprintf(fdata, "# channel estimation: cunning\n");
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

	fprintf(fdata, "# Relay EbN0: Fixed to %d dB \n", RLY_CNR);
	printf("# Relay EbN0: Fixed to %d dB \n", RLY_CNR);

#if(RS_DF==ON)
	fprintf(fdata, "# Relay Signalling: DF\n");
	printf("# Relay Signalling: DF\n");
#else
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
#if((SWITCH==ON)&&(RS_DF==OFF))
	cout << "Error! switch on and its AF" << endl;
#elif((SWITCH==OFF)&&(RS_DF==ON))
	cout << "Error! switch off and its DF" << endl;
#endif

#ifdef  UD_SELECT	
#ifndef	E_SELECT
	fprintf(fdata, "# IWI selection: ON\n");
	printf("# IWI selection: ON\n");
#else							  
	fprintf(fdata, "# IEI selection: ON\n");
	printf("# IEI selection: ON\n");
#endif
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
	/*
	fprintf(fdata, "# Error: Use HD for MRC\n");
	printf("# Error: Use HD for MRC\n");
	_getch();  */
	fprintf(fdata, "# Combining: MRC\n");
	printf("# Combining: MRC\n");
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
	fprintf(fdata, "# CD: ERROR!(SDVD)\n");
	printf("# CD: SDVD\n");
#endif

#ifdef HD
	fprintf(fdata, "# CD: HDVD\n");
	printf("# CD: HDVD\n");
#endif						
	fprintf(fdata, "# CELL EDGE: %d dB \n", LOSS);
	printf("# CELL EDGE: %d dB \n", LOSS); 
#if(RLY_EXPMODEL==ON)
	cout << "Exponential Channel Profile" << endl;
	cout << "Tau_rs: " << RS_TAU << endl;
	cout << "Tau_ms: " << TAU << endl;	 
#endif					

	fclose(fdata);			
}
