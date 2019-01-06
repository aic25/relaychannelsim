#include "const.h"
double h_r[BS][PATH][2];
void  relay(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs, double(*BS_signal_power)){

#if SWITCH == 0
	int k;
	for (k = 0; k < BURST; k++){
		relay_out[k][0] = relay_in[k][0];
		relay_out[k][1] = relay_in[k][1];
	}

#elif SWITCH == 1
#ifdef DF
	Relay_DF(relay_in, relay_out, bit, bs);
#else
	Relay_AF(relay_in, relay_out, bit, bs, BS_signal_power);
#endif // DF	 
#endif
}

void Relay_DF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs){
	double relay_received_signal[BURST][2];
	double RS_noise_power[BS];
	int k;

	RS_noise_power[bs] = 0;

#ifdef RLY_PROPOGATION
	Relay_propagation(relay_in, relay_received_signal, bs);	/* Signal goes thru channel and reaches the relay element */
#else
	for ( k = 0; k < BURST; k++){
		relay_received_signal[k][0] = relay_in[k][0];
		relay_received_signal[k][1] = relay_in[k][1];
	}
#endif // RLY_PROPOGATION

#ifdef RLY_NOISE
#ifdef RLY_PROPOGATION
	noise_addition(0, relay_received_signal, relay_received_signal, &RS_noise_power[bs]);	/* given the CNR add noise */
#else
	noise_addition(RLY_CNR, relay_received_signal, relay_received_signal, &RS_noise_power[]);
#endif // RLY_PROPOGATION
#endif // RLY_NOISE

	Relay_receiver(relay_received_signal, bit, bs);	/* HD decoder. */
	Relay_transmitter(bit, relay_out);	/* Encode, modulate etc. */
}

void Relay_AF(double(*relay_in)[2], double(*relay_out)[2], int(*bit)[N_dbps], int bs, double(*BS_signal_power)){
	double RS_signal_power[BS], beta;
	double RS_noise_power[BS];
	int k;

	RS_noise_power[bs] = 0;
	RS_signal_power[bs] = 0;

#ifdef RLY_PROPOGATION
	Relay_propagation(relay_in, relay_out, bs);	/* Signal goes thru channel and reaches the relay element */
#else
	for ( k = 0; k < BURST; k++){
		relay_out[k][0] = relay_in[k][0];
		relay_out[k][1] = relay_in[k][1];
	}
#endif // RLY_PROPOGATION

#ifdef RLY_NOISE
#ifdef RLY_PROPOGATION
	noise_addition(0, relay_out, relay_out, &RS_noise_power[bs]);	/* given the CNR add noise */
#else
	noise_addition(RLY_CNR, relay_out, relay_out, &RS_noise_power[bs]);
#endif // RLY_PROPOGATION
#endif // RLY_NOISE				

	/* Beta is calculated */
	for ( k = 0; k < BURST; k++){
		RS_signal_power[bs] += sqr(relay_out[k][0]) + sqr(relay_out[k][1]);
	}
	RS_signal_power[bs] /= BURST;

	beta = sqrt(BS_signal_power[bs] / RS_signal_power[bs]);

	for ( k = 0; k < BURST; k++){
		relay_out[k][0] = beta*relay_out[k][0];
		relay_out[k][1] = beta*relay_out[k][1];
	}
}
/************************* Relay Propagation and Noise addition *******************************/

void Relay_propagation(double(*transmitted_signal)[2], double(*Received_Signal)[2], int bs){
	int time, d, index;

	for (time = 0; time < BURST; time++) {
		Relay_fading_process(time, h_r[bs]);	/* Channel generator */

		Received_Signal[time][0] = 0.0;
		Received_Signal[time][1] = 0.0;

		for (d = 0; d < PATH; d++) {
			index = time - d * DELAY;
			if (0 <= index && index < BURST) {
				Received_Signal[time][0] += h_r[bs][d][0] * transmitted_signal[index][0] - h_r[bs][d][1] * transmitted_signal[index][1];
				Received_Signal[time][1] += h_r[bs][d][1] * transmitted_signal[index][0] + h_r[bs][d][0] * transmitted_signal[index][1];
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





