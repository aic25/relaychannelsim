#include "const.h"		  
extern double CNR, B2_signal_power;
double h_r2[PATH][2];

void  relay2(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]){

#if SWITCH == 0
	int k;
	for(k=0; k<BURST; k++){
		relay_out[k][0] = relay_in[k][0];
		relay_out[k][1] = relay_in[k][1];
	}	  

#elif SWITCH == 1
#ifdef DF
	Relay2_DF(relay_in, relay_out, bit);
#else
	Relay2_AF(relay_in, relay_out, bit);
#endif // DF	
#endif
}

void  Relay2_DF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]){
	double relay_received_signal[BURST][2];

#ifdef RLY_PROPOGATION
	Relay2_propagation(relay_in, relay_received_signal);
#else
	for (size_t k = 0; k < BURST; k++){
		relay_received_signal[k][0] = relay_in[k][0];
		relay_received_signal[k][1] = relay_in[k][1];
	}
#endif // RLY_PROPOGATION

#ifdef RLY_NOISE
#ifdef RLY_PROPOGATION
	noise_addition(0, relay_received_signal, relay_received_signal);	/* given the CNR add noise */
#else
	noise_addition(RLY_CNR, relay_received_signal, relay_received_signal);
#endif // RLY_PROPOGATION
#endif // RLY_NOISE

	Relay2_receiver(relay_received_signal, bit);	/* HD decoder. */
	Relay2_transmitter(bit, relay_out);	/* Encode, modulate etc. */
}

void Relay2_AF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]){
	double R2_signal_power, beta;

	R2_signal_power = 0;

#ifdef RLY_PROPOGATION
	Relay2_propagation(relay_in, relay_out);	/* Signal goes thru channel and reaches the relay element */
#else
	for (size_t k = 0; k < BURST; k++){
		relay_out[k][0] = relay_in[k][0];
		relay_out[k][1] = relay_in[k][1];
	}
#endif // RLY_PROPOGATION

#ifdef RLY_NOISE
#ifdef RLY_PROPOGATION
	noise_addition(0, relay_out, relay_out);	/* given the CNR add noise */
#else
	noise_addition(RLY_CNR, relay_out, relay_out);
#endif // RLY_PROPOGATION
#endif // RLY_NOISE				

	/* Beta is calculated */
	for (size_t k = 0; k < BURST; k++){
		R2_signal_power += sqr(relay_out[k][0]) + sqr(relay_out[k][1]);
	}
	R2_signal_power /= BURST;

	beta = sqrt(B2_signal_power / R2_signal_power);

	for (size_t k = 0; k < BURST; k++){
		relay_out[k][0] = beta*relay_out[k][0];
		relay_out[k][1] = beta*relay_out[k][1];
	}
}

void Relay2_propagation(double(*transmitted_signal)[2], double(*received_signal)[2]){
	int time, d, index;

	for (time = 0; time < BURST; time++) {
		Relay_fading_process(time, h_r2);

		received_signal[time][0] = 0.0;
		received_signal[time][1] = 0.0;

		for (d = 0; d < PATH; d++) {
			index = time - d * DELAY;
			if (0 <= index && index < BURST) {
				received_signal[time][0] += h_r2[d][0] * transmitted_signal[index][0] - h_r2[d][1] * transmitted_signal[index][1];
				received_signal[time][1] += h_r2[d][1] * transmitted_signal[index][0] + h_r2[d][0] * transmitted_signal[index][1];
			}
		}
	}
}
