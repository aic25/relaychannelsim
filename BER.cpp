#include "const.h"	
extern double RLY_CNR, CNR, EbNo, signal_power, noise_power;
extern int maxloop;
extern errno_t err;

void BER_RsPlusBs(int loop, int(*rs_bit)[N_dbps], int(*tbit)[N_dbps], int(*rbit)[N_dbps]){
	FILE *fdata, *output_log, *rs_fdata;
	std::string file_time;
	int	i, k, bitError = 0, rs_bitError = 0;
	double Pe, rs_Pe;
	static double AverageBER = 0.0, rs_AverageBER = 0.0;
	static int pktError = 0, rs_pktError = 0;


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
			if (tbit[k][i] != rs_bit[k][i]) rs_bitError++;
		}
	}

	rs_Pe = (double)rs_bitError / (Nd*(N_dbps - TAIL));
	Pe = (double)bitError / (Nd*(N_dbps - TAIL));

	rs_AverageBER += (double)rs_Pe / (maxloop);
	AverageBER += (double)Pe / (maxloop);

	if (rs_bitError != 0) rs_pktError++;
	if (bitError != 0) pktError++;

	printf("Eb/N0=%2d::loop=%4d:rs_BER = %e :BER = %e %e \n", (int)CNR - 3 * MODULATION / 2, loop, rs_Pe, Pe, AverageBER);
	fprintf(output_log, "Eb/N0=%2d, loop=%3d, rs_BER = %e, BER = %e %e \n", (int)CNR - 3 * MODULATION / 2, loop, rs_Pe, Pe, AverageBER);

	if (loop == maxloop - 1){
		printf("Eb/No = %ld(dB) :: %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)CNR - 3 * MODULATION / 2, 10.0*log10(signal_power / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));
		fprintf(output_log, "Eb/No = %ld(dB), %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)CNR - 3 * MODULATION / 2, 10.0*log10(signal_power / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));

		if ((err = fopen_s(&fdata, "BER.dat", "a")) != 0)	{
			printf("The file BER.dat file was not opened\n");
			getchar();
			exit(1);
		};
		if ((err = fopen_s(&rs_fdata, "rs_BER.dat", "a")) != 0)	{
			printf("The file rs_BER.dat file was not opened\n");
			getchar();
			exit(1);
		};
		fprintf(fdata, "%1d\t%e\n", (int)CNR - 3 * MODULATION / 2, AverageBER);
		fprintf(rs_fdata, "%e\n", rs_AverageBER);
		fclose(fdata);
		fclose(rs_fdata);
		AverageBER = 0.0;
		pktError = 0;
	}
	fclose(output_log);
}

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
	printf("Eb/N0=%2d :: loop=%3d :: BER = %e %e \n", (int)RLY_CNR - 3 * MODULATION / 2, loop, Pe, AverageBER);
	fprintf(output_log, "Eb/N0=%2d, loop=%3d, BER = %e %e \n", (int)RLY_CNR - 3 * MODULATION / 2, loop, Pe, AverageBER);

	if (loop == maxloop - 1){
		printf("Eb/No = %ld(dB) :: %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)RLY_CNR - 3 * MODULATION / 2, 10.0*log10(signal_power / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));
		fprintf(output_log, "Eb/No = %ld(dB), %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)RLY_CNR - 3 * MODULATION / 2, 10.0*log10(signal_power / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));

		if ((err = fopen_s(&fdata, "BER.dat", "a")) != 0)	{
			printf("The file BER.dat file was not opened\n");
			getchar();
			exit(1);
		};
		fprintf(fdata, "%1d\t%e\n", (int)RLY_CNR - 3 * MODULATION / 2, AverageBER);
		fclose(fdata);
		AverageBER = 0.0;
		pktError = 0;
	}
	fclose(output_log);
}