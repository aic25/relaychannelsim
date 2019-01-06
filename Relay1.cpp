#include "const.h"
extern double CNR, B1_signal_power;
double h_r[PATH][2];
void  relay1(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]){

#if SWITCH == 0
	int k;
	for (k = 0; k < BURST; k++){
		relay_out[k][0] = relay_in[k][0];
		relay_out[k][1] = relay_in[k][1];
	}

#elif SWITCH == 1
#ifdef DF
	Relay1_DF(relay_in, relay_out, bit);
#else
	Relay1_AF(relay_in, relay_out, bit);
#endif // DF	 
#endif
}

void Relay1_DF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]){
	double relay_received_signal[BURST][2];

#ifdef RLY_PROPOGATION
	Relay1_propagation(relay_in, relay_received_signal);	/* Signal goes thru channel and reaches the relay element */
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

	Relay1_receiver(relay_received_signal, bit);	/* HD decoder. */
	Relay1_transmitter(bit, relay_out);	/* Encode, modulate etc. */
}

void Relay1_AF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps]){
	double R1_signal_power, beta;

	R1_signal_power = 0;

#ifdef RLY_PROPOGATION
	Relay1_propagation(relay_in, relay_out);	/* Signal goes thru channel and reaches the relay element */
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
		R1_signal_power += sqr(relay_out[k][0]) + sqr(relay_out[k][1]);
	}
	R1_signal_power /= BURST;

	beta = sqrt(B1_signal_power / R1_signal_power);

	for (size_t k = 0; k < BURST; k++){
		relay_out[k][0] = beta*relay_out[k][0];
		relay_out[k][1] = beta*relay_out[k][1];
	}
}
/************************* Relay Propagation and Noise addition *******************************/

void Relay1_propagation(double(*transmitted_signal)[2], double(*received_signal)[2]){
	int time, d, index;

	for (time = 0; time < BURST; time++) {
		Relay_fading_process(time, h_r);	/* Channel generator */

		received_signal[time][0] = 0.0;
		received_signal[time][1] = 0.0;

		for (d = 0; d < PATH; d++) {
			index = time - d * DELAY;
			if (0 <= index && index < BURST) {
				received_signal[time][0] += h_r[d][0] * transmitted_signal[index][0] - h_r[d][1] * transmitted_signal[index][1];
				received_signal[time][1] += h_r[d][1] * transmitted_signal[index][0] + h_r[d][0] * transmitted_signal[index][1];
			}
		}
	}
}

void Relay_fading_process(int time, double(*w)[2]){
	int	d, kk;
	double ftemp, phase, x, power[PATH][RAY_COMPONENT];
	static double omega[PATH][RAY_COMPONENT], p0[PATH][RAY_COMPONENT], a[PATH][RAY_COMPONENT];

	if (time == 0) {
		ftemp = 2.0 * PI / SAMPLING_RATE;
		for (d = 0; d < PATH; d++){
			for (kk = 0; kk < RAY_COMPONENT; kk++) {
				omega[d][kk] = ftemp * FDr * cos(2.0 * PI * (double)rand() / RAND_MAX);
				p0[d][kk] = 2.0 * PI * (double)rand() / RAND_MAX;
				if (EXPMODEL == ON) {
					exponential_PATH(d, kk, power);
				}
				else{
					samelevel_PATH(d, kk, power);
				}
				x = (double)rand() / RAND_MAX;
				if (x < 1.0e-6) x = 1.0e-6;
				a[d][kk] = sqrt(power[d][kk] * -log(x) * (pow(10.0, RLY_CNR / 10.0) / CODING_RATE));
			}
		}
	}
	for (d = 0; d < PATH; d++) {
		w[d][0] = 0.0;
		w[d][1] = 0.0;

		for (kk = 0; kk < RAY_COMPONENT; kk++) {
			phase = omega[d][kk] * (double)time + p0[d][kk];

			w[d][0] += a[d][kk] * cos(phase);
			w[d][1] += a[d][kk] * sin(phase);
		}
	}
}





