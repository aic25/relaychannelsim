#include "const-function-declerations.h"	
extern double CNR, EbNo, signal_power, noise_power;
extern int maxloop, R_index[N];
extern errno_t err;

void BER_RsPlusBs(int loop, int(*rs_bit)[Nd][N_dbps], int(*tbit)[Nd][N_dbps], int(*rbit)[Nd][N_dbps]){
	FILE *fdata, *output_log, *rs_fdata;
	std::string file_time;
	int	i, k, bitError = 0, rs_bitError = 0;
	double Pe, rs_Pe;
	static double AverageBER = 0.0, rs_AverageBER = 0.0;
	static int pktError = 0, rs_pktError = 0;


	if ((err = fopen_s(&output_log, "output_log.dat", "a")) != 0)	{ printf("The file outputlog file was not opened\n");	getchar();	exit(1); };

	file_time = currentDateTime();
	fprintf(output_log, "# %s ==  ", file_time.c_str());

	for (k = 0; k < Nd; k++){
		for (i = 0; i < N_dbps - TAIL; i++) {
			if ((tbit[0][k][i] != rbit[0][k][i])) bitError++;
			if ((RS_DF == ON) && ((tbit[0][k][i] != rs_bit[0][k][i]) || (tbit[1][k][i] != rs_bit[1][k][i]))) rs_bitError++;
		}
	}

	rs_Pe = (double)rs_bitError / (KI*Nd*(N_dbps - TAIL));
	Pe = (double)bitError / (Nd*(N_dbps - TAIL));

	rs_AverageBER += (double)rs_Pe / (maxloop);
	AverageBER += (double)Pe / (maxloop);

	if (rs_bitError != 0) rs_pktError++;
	if (bitError != 0) pktError++;

	printf("Eb/N0=%2d::loop=%5d:rs_BER = %1.2e :BER = %1.2e %1.2e \n", (int)CNR - 3 * MODULATION / 2, loop, rs_Pe, Pe, AverageBER);
	fprintf(output_log, "Eb/N0=%2d, loop=%5d, rs_BER = %1.2e, BER = %1.2e %1.2e \n", (int)CNR - 3 * MODULATION / 2, loop, rs_Pe, Pe, AverageBER);

	if (loop == maxloop - 1){
		printf("Eb/No = %f(dB) :: %f(dB)\tPER = %e\tAveBER = %e\n",
			CNR - 3 * MODULATION / 2, 10.0*log10((signal_power) / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));
		fprintf(output_log, "Eb/No = %f(dB), %f(dB)\tPER = %e\tAveBER = %e\n",
			CNR - 3 * MODULATION / 2, 10.0*log10((signal_power) / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));

		if ((err = fopen_s(&fdata, "BER.dat", "a")) != 0)	{ printf("The file BER.dat file was not opened\n");	getchar();	exit(1); };
		if ((err = fopen_s(&rs_fdata, "rs_BER.dat", "a")) != 0)	{ printf("The file rs_BER.dat file was not opened\n");	getchar();	exit(1); };

		fprintf(fdata, "%f\t%e\n", CNR - 3 * MODULATION / 2, AverageBER);
		fprintf(rs_fdata, "%e\n", rs_AverageBER);
		fclose(fdata);
		fclose(rs_fdata);
		AverageBER = 0.0;
		pktError = 0;
	}
	fclose(output_log);
}

