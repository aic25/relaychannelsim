#include "const.h"	
extern errno_t err;

void BER(int loop, int maxloop, double CNR, int(*tbit)[Nd][N_dbps], int(*rbit)[Nd][N_dbps], double noise_power, double signal_power, double interference_power){
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
			if ((tbit[0][k][i] != rbit[0][k][i]) || (tbit[1][k][i] != rbit[1][k][i])) bitError++;
		}
	}

	Pe = (double)bitError / (2*Nd*(N_dbps - TAIL));
	AverageBER += (double)Pe / (maxloop);
	if (bitError != 0) pktError++;
	printf("Eb/N0=%2d :: loop=%3d :: BER = %e %e \n", (int)CNR - 3 * MODULATION / 2, loop, Pe, AverageBER);
	fprintf(output_log, "Eb/N0=%2d, loop=%3d, BER = %e %e \n", (int)CNR - 3 * MODULATION / 2, loop, Pe, AverageBER);
	
	if (loop == maxloop - 1){
		printf("Eb/No = %ld(dB) :: %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)CNR - 3 * MODULATION / 2, 10.0*log10(signal_power / (CODING_RATE*noise_power)),			// -3 deleted
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));
		fprintf(output_log, "Eb/No = %ld(dB), %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)CNR - 3 * MODULATION / 2, 10.0*log10(signal_power / (CODING_RATE*noise_power)),
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));	 

		if ((err = fopen_s(&fdata, "BER.dat", "a")) != 0)	{
			printf("The file BER.dat file was not opened\n");
			getchar();
			exit(1);
		};	   
		fprintf(fdata, "%1d\t%e\n", (int)CNR - 3 * MODULATION / 2, AverageBER);
		fclose(fdata);
		AverageBER = 0.0;
		pktError = 0;
	}
	fclose(output_log);
}	