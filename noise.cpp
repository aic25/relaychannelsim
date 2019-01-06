#include "const.h"
extern double noise_power;
void noise_addition(double CNR, double(*input_signal)[2], double(*output_signal)[2]){
	int	time_index;
	double noise[2], noise_amplitude;

	noise_amplitude = sqrt(pow(10.0, -CNR / 10.0) / CODING_RATE);

	for (time_index = 0; time_index < BURST; time_index++) {
		noise_generator(noise_amplitude, noise);
		output_signal[time_index][0] = input_signal[time_index][0] + noise[0];
		output_signal[time_index][1] = input_signal[time_index][1] + noise[1];
	}
}

void noise_generator(double n_amplitude, double(*noise)){
	double a, b, sqb;

	a = (double)rand() / RAND_MAX;
	b = (double)rand() / RAND_MAX;

	if (b < 1.0e-6) b = 1.0e-6;
	sqb = sqrt(-log(b)) * n_amplitude;

	noise[0] = sqb * cos(2.0 * PI * a);
	noise[1] = sqb * sin(2.0 * PI * a);
	noise_power += sqr(noise[0]) + sqr(noise[1]);
}

