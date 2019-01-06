#include "const.h"	
extern double CNR, EbNo, signal_power, noise_power;
extern int XI;
extern int maxloop;
extern errno_t err;

void BER(int loop, int(*tbit)[N_dbps], int(*rbit)[N_dbps]){
	FILE *fdata, *output_log;
	std::string file_time;
	int	i, k, bitError = 0;
	double Pe;
	static double AverageBER = 0.0;
	static int pktError = 0;


	if ((err = fopen_s(&output_log, "output_log.dat", "a")) != 0)	{
		printf("The file outputlog file was not opened\n");
		getchar();
		exit(1);
	};

	file_time = currentDateTime();
	fprintf(output_log, "# %s ==  ", file_time.c_str());

	for (k = 0; k < Nd; k++){
		for (i = 0; i < N_dbps - TAIL; i++) {
			if (tbit[k][i] != rbit[k][i]) bitError++;
		}
	}

	Pe = (double)bitError / (Nd*(N_dbps - TAIL));
	AverageBER += (double)Pe / (maxloop);
	if (bitError != 0) pktError++;
	printf("XI=%2d :: loop=%3d :: BER = %e %e \n", (int)XI, loop, Pe, AverageBER);
	fprintf(output_log, "XI=%2d, loop=%3d, BER = %e %e \n", (int)XI, loop, Pe, AverageBER);
	
	if (loop == maxloop - 1){
		printf("XI = %ld(dB) :: %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)XI, 10.0*log10(signal_power / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));
		fprintf(output_log, "XI = %ld(dB), %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)XI, 10.0*log10(signal_power / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));	 

		if ((err = fopen_s(&fdata, "BER.dat", "a")) != 0)	{
			printf("The file BER.dat file was not opened\n");
			getchar();
			exit(1);
		};	   
		fprintf(fdata, "%1d\t%e\n", (int)XI, AverageBER);
		fclose(fdata);
		AverageBER = 0.0;
		pktError = 0;
	}
	fclose(output_log);
}

