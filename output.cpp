#include "const.h"	
extern errno_t err;
extern double Pe;

void BER(int(*R_index), int loop, int maxloop, double CNR, int(*tbit)[Nd][N_dbps], int(*rbit)[Nd][N_dbps], double noise_power, double signal_power, double interference_power){
	FILE *fdata, *output_log;
	std::string file_time;
	int	i, k, bitError = 0;
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
			if ((INTERFERENCE == ON) ? ((tbit[0][k][i] != rbit[0][k][i]) || (tbit[1][k][i] != rbit[1][k][i])) : (tbit[0][k][i] != rbit[0][k][i])) bitError++;
		}
	}

	Pe = (double)bitError / ((INTERFERENCE == ON) ? (2 * Nd*(N_dbps - TAIL)) : (Nd*(N_dbps - TAIL)));
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

void Spatio_Temporal_Printout(int loop, double(*avg), double(*e)[NUMBER_OF_STEPS][2], char File_Name[15]){
	FILE *F_pointer;
	char strI[10];
	char Ext[5] = ".dat";
	int step, f_index;

	_itoa(loop, strI, 10);
	strcat(File_Name, strI);
	strcat(File_Name, Ext);
	if ((err = fopen_s(&F_pointer, File_Name, "w+")) != 0)	{
		printf("The file F_pointer.dat file was not opened\n");
		getchar();
		exit(1);
	};
	//File_Name[4] = '\0';
	//printf("%d\t%d\n",loop,EbNo);
	for (step = 0; step < Np + Nd; step++)
	{
		avg[step] = 0;
		for (f_index = 0; f_index < N; f_index++)
		{
			avg[step] += sqrt(sqr(e[f_index][step][0]) + sqr(e[f_index][step][1]));
			fprintf(F_pointer, "%d", step);
			fprintf(F_pointer, "\t%f", sqrt(sqr(e[f_index][step][0]) + sqr(e[f_index][step][1])));
		}
		avg[step] = avg[step] / N;
		//fprintf(F_pointer, "%d\t%e\n", step, e_rls_avg[step]);
		fprintf(F_pointer, "\n");
	}
	fclose(F_pointer);
}

double Take_Average(double(*list), int Number_of_Elements, int From){
	double sum;
	int i;
	sum = 0;
	for (i = From; i < Number_of_Elements; i++){
		sum += list[i];
	}
	sum /= Number_of_Elements;
	return sum;
}

double Take_Average_C(double(*list)[2], int Number_of_Elements, int From){
	double sum;
	int i;
	sum = 0;
	for (i = From; i < Number_of_Elements; i++){
		sum += sqrt(sqr(list[i][0]) + sqr(list[i][1]));
	}
	sum /= Number_of_Elements;
	return sum;
}

void Correlation(double EbN0, int loop, double(*e)[NUMBER_OF_STEPS][2], double(*Criterion)[NUMBER_OF_STEPS][2], char File_Name[25]){
	FILE *F_pointer;
	double e_av[N], C_av[N], tmp, C_av_av, e_av_av;
	int i, step;
	for (i = 0; i < N; i++){
		e_av[i] = Take_Average_C(e[i], (int)NUMBER_OF_STEPS, 8);
	}
	for (i = 0; i < N; i++){
		C_av[i] = Take_Average_C(Criterion[i], (int)NUMBER_OF_STEPS, 8);
	}

	C_av_av = Take_Average(C_av, (int)N, 0);
	e_av_av = Take_Average(e_av, (int)N, 0);

	for (i = 0; i < N; i++){
		C_av[i] -= C_av_av;
		e_av[i] -= e_av_av;
	}

	char Ext[5] = ".dat";
	char ebn0[10];

#ifdef SEPARATE_CORR_PLOOP 
	char strI[10];
	_itoa(loop, strI, 10);
	strcat(File_Name, strI);	
	strcat(File_Name, "_");
#endif // DEBUG

	_itoa(EbN0, ebn0, 10);
	strcat(File_Name, ebn0);
	strcat(File_Name, Ext);
	if ((err = fopen_s(&F_pointer, File_Name, "a+")) != 0)	{
		printf("The file F_pointer.dat file was not opened\n");
		getchar();
		exit(1);
	};

	for (step = 0; step < N; step++){
		//fprintf(F_pointer, "%d\t%f\t%f\n", step, C_av[step], e_av[step]);
		fprintf(F_pointer, "%f\t%f\n", e_av[step], C_av[step]);
	}
	fclose(F_pointer);

}