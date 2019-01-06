#include "const-function-declerations.h"	
extern double CNR, EbNo, signal_power, noise_power;
extern int maxloop, R_index[N];
extern errno_t err;
double theoretical_BER;

void BER_RsPlusBs(int loop, int(*rs_bit)[Nd][N_dbps], int(*tbit)[Nd][N_dbps], int(*rbit)[Nd][N_dbps]) {
	FILE *fdata, *output_log, *rs_fdata;
	std::string file_time;
	int	i, k, bitError = 0, rs_bitError = 0;
	double Pe, rs_Pe;
	static double AverageBER = 0.0, rs_AverageBER = 0.0;
	static int pktError = 0, rs_pktError = 0;


	if ((err = fopen_s(&output_log, "output_log.dat", "a")) != 0) { printf("The file outputlog file was not opened\n");	getchar();	exit(1); };

	file_time = currentDateTime();
	fprintf(output_log, "# %s ==  ", file_time.c_str());

	for (k = 0; k < Nd; k++) {
		for (i = 0; i < N_dbps - TAIL; i++) {
			if ((tbit[0][k][i] != rbit[0][k][i]) || (tbit[1][k][i] != rbit[1][k][i])) bitError++;
			if ((RS_DF == ON) && ((tbit[0][k][i] != rs_bit[0][k][i]) || (tbit[1][k][i] != rs_bit[1][k][i]))) rs_bitError++;
		}
	}

	rs_Pe = (double)rs_bitError / (KI*Nd*(N_dbps - TAIL));
	Pe = (double)bitError / (KI*Nd*(N_dbps - TAIL));

	rs_AverageBER += (double)rs_Pe / (maxloop);
	AverageBER += (double)Pe / (maxloop);

	if (rs_bitError != 0) rs_pktError++;
	if (bitError != 0) pktError++;

	printf("Eb/N0=%2d::loop=%5d:rs_BER = %1.2e :BER = %1.2e %1.2e \n", (int)CNR - 3 * MODULATION / 2, loop, rs_Pe, Pe, AverageBER);
	fprintf(output_log, "Eb/N0=%2d, loop=%5d, rs_BER = %1.2e, BER = %1.2e %1.2e \n", (int)CNR - 3 * MODULATION / 2, loop, rs_Pe, Pe, AverageBER);

	if (loop == maxloop - 1) {
		printf("Eb/No = %ld(dB) :: %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)CNR - 3 * MODULATION / 2, 10.0*log10((signal_power) / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));
		fprintf(output_log, "Eb/No = %ld(dB), %f(dB)\tPER = %e\tAveBER = %e\n",
			(int)CNR - 3 * MODULATION / 2, 10.0*log10((signal_power) / (CODING_RATE*noise_power)) - 3.0,
			(double)pktError / (double)(loop + 1),
			AverageBER / (double)(loop + 1));

		if ((err = fopen_s(&fdata, "BER.dat", "a")) != 0) { printf("The file BER.dat file was not opened\n");	getchar();	exit(1); };
		if ((err = fopen_s(&rs_fdata, "rs_BER.dat", "a")) != 0) { printf("The file rs_BER.dat file was not opened\n");	getchar();	exit(1); };

		fprintf(fdata, "%1d\t%e\n", (int)CNR - 3 * MODULATION / 2, AverageBER);
		fprintf(rs_fdata, "%e\n", rs_AverageBER);
		fclose(fdata);
		fclose(rs_fdata);
		AverageBER = 0.0;
		pktError = 0;
	}
	fclose(output_log);
}

void Theoretical_BER(int loop, double(*H)[TS][BS][2], double(*W)[BS][CL][2]) {
	FILE *theoretical_ber_data;
	std::string file_time;
	int f_index, sb1hat, sb1, sb2hat, sb2; /* counters */
	double Nb = MODULATION;
	double
		Ne,
		Pe,
		sigma_nu2,
		sigma_nu,
		x, /* nominator */
		x_max,
		x_min,
		denominator, /* denominator of erfc */
		theo_BER,
		theo_BER_N,
		sigma_rk2 = pow(10.0, -5),	  /* sigma_rk is assumed to be the same for all sb1 */
		sigma_m121 = pow(10.0, -CNR / 10.0),
		sigma_m122 = pow(10.0, -CNR / 10.0); /* LOSS is applied to the channel coefficient therefore noise power is the same between phases */
	std::complex < double >
		temp,
		a_1,
		a_2,
		a_3,
		a_4,
		w_0c,
		w_1c,
		w_2c,
		h_m1b1,
		h_m1b2,
		h_m1b3,
		h_m1r1,
		h_m1r2,
		h_m1r3,
		DeltaSb1,
		DeltaSb2;
	/*int Delta_Sb_1[16][3] = {
		{ 0,0,0 },{ 0,-2,1 },{ -2,-2,2 },{ -2,0,1 },
		{ 0,2,1 },{ 0,0,0 },{ -2,0,1 },{ -2,2,2 },
		{ 2,2,2 },{ 2,0,1},{ 0,0,0 },{ 0,2,1 },
		{ 2,0,1 },{ 2,-2,2 },{ 0,-2,1 },{ 0,0,0 }
		}; /* third collumn number of erroneous bits */
	/*int Delta_Sb_2[16][2] = {
		{0,0},{0,-2},{-2,-2},{-2,0},
		{0,2},{0,0},{-2,0},{-2,2},
		{2,2},{2,0},{0,0},{0,2},
		{2,0},{2,-2},{0,-2},{0,0}
		}; /* \Delta{S}_{b_2} = */
	Vector4cd Symbols = {
		complex<double>(-OneBySqrt2, -OneBySqrt2),
		complex<double>(-OneBySqrt2, OneBySqrt2),
		complex<double>(OneBySqrt2, OneBySqrt2),
		complex<double>(OneBySqrt2, -OneBySqrt2),
	};
	/*cout << "(00)-(01): " << abs(complex<double>(-OneBySqrt2, -OneBySqrt2) - complex<double>(-OneBySqrt2, OneBySqrt2)) << endl;
	cout << "(00)-(01): " << abs(complex<double>(-OneBySqrt2, -OneBySqrt2) - complex<double>(OneBySqrt2, OneBySqrt2)) << endl;
	cout << "E{S^2}: " << pow(complex<double>(-OneBySqrt2, -OneBySqrt2), 2) << endl;
	cout << "E{|S|^2}: " << pow(abs(complex<double>(-OneBySqrt2, -OneBySqrt2)), 2) << endl;
	cout << "S*S^*: " << complex<double>(-OneBySqrt2, -OneBySqrt2)*complex<double>(-OneBySqrt2, OneBySqrt2) << endl;
	cout << "S*S: " << complex<double>(-OneBySqrt2, -OneBySqrt2)*complex<double>(-OneBySqrt2, -OneBySqrt2) << endl;
	getchar();*/

	loop == 0 ? theoretical_BER = 0 : theoretical_BER = theoretical_BER;
	theo_BER_N = 0;
	for (f_index = 0; f_index < N; f_index++)
		/* Loop length is over all subarriers.
		After the average over all subcarriers is obtained
		an average value should be updated with each value untill maxloop. */
	{
		h_m1b1 = complex<double>(H[f_index][0][0][0], H[f_index][0][0][1]);
		h_m1b2 = complex<double>(H[f_index][0][1][0], H[f_index][0][1][1]);
		h_m1b3 = complex<double>(H[f_index][0][2][0], H[f_index][0][2][1]);
		h_m1r1 = complex<double>(H[f_index][1][0][0], H[f_index][1][0][1]);
		h_m1r2 = complex<double>(H[f_index][1][1][0], H[f_index][1][1][1]);
		h_m1r3 = complex<double>(H[f_index][1][2][0], H[f_index][1][2][1]);
		w_0c = complex<double>(W[f_index][0][0][0], -W[f_index][0][0][1]);
		w_1c = complex<double>(W[f_index][0][1][0], -W[f_index][0][1][1]);
		w_2c = complex<double>(W[f_index][0][2][0], -W[f_index][0][2][1]);
		/*calculation of sigma_nu */
		sigma_nu2
			= pow(abs(w_1c*h_m1b3 + w_2c*h_m1r3), 2)
			+ pow(abs(w_1c), 2)*sigma_m121 + pow(abs(w_2c), 2)*sigma_m122
			+ pow(abs(w_2c*h_m1r1), 2)*sigma_rk2
			+ pow(abs(w_2c*h_m1r2), 2)*sigma_rk2
			+ pow(abs(w_2c*h_m1r3), 2)*sigma_rk2;/**/
		sigma_nu = sqrt(sigma_nu2);
		/* calculation of a_x */
		a_1 = 1;
		a_2 = w_0c;
		a_3 = -(w_1c*h_m1b1 + w_2c*h_m1r1);
		a_4 = -(w_1c*h_m1b2 + w_2c*h_m1r2);
		/* cout << "W1:" << endl;
		cout << W[f_index][0][0][0] << "," << W[f_index][0][0][1] << endl;
		cout << W[f_index][0][1][0] << "," << W[f_index][0][1][1] << endl;
		cout << W[f_index][0][2][0] << "," << W[f_index][0][2][1] << endl;
		cout << "W2:" << endl;
		cout << W[f_index][1][0][0] << "," << W[f_index][1][0][1] << endl;
		cout << W[f_index][1][1][0] << "," << W[f_index][1][1][1] << endl;
		cout << W[f_index][1][2][0] << "," << W[f_index][1][2][1] << endl;
		cout << "W3:" << endl;
		cout << W[f_index][2][0][0] << "," << W[f_index][2][0][1] << endl;
		cout << W[f_index][2][1][0] << "," << W[f_index][2][1][1] << endl;
		cout << W[f_index][2][2][0] << "," << W[f_index][2][2][1] << endl;*/
		/*cout << "H" << endl
			<< complex<double>(H[f_index][0][0][0], H[f_index][0][0][1]) << ", "
			<< complex<double>(H[f_index][1][0][0], H[f_index][1][0][1]) << endl
			<< complex<double>(H[f_index][0][1][0], H[f_index][0][1][1]) << ", "
			<< complex<double>(H[f_index][1][1][0], H[f_index][1][1][1]) << endl
			<< complex<double>(H[f_index][0][2][0], H[f_index][0][2][1]) << ", "
			<< complex<double>(H[f_index][1][2][0], H[f_index][1][2][1]) << endl;
			cout << "Sigma_nu2 part 1: "
			<< pow(abs(complex<double>(W[f_index][0][0][0], -W[f_index][0][0][1])
			*complex<double>(H[f_index][0][2][0], H[f_index][0][2][1])
			+ complex<double>(W[f_index][0][1][0], -W[f_index][0][1][1])
			*complex<double>(H[f_index][1][2][0], H[f_index][1][2][1])), 2) << endl;
			cout << "Sigma_nu2 part 2: "
			<< +pow(abs(complex<double>(W[f_index][0][0][0], -W[f_index][0][0][1])), 2)*sigma_m121
			+ pow(abs(complex<double>(W[f_index][0][1][0], -W[f_index][0][1][1])), 2)*sigma_m122 << endl;
			cout << "Sigma_nu2 part 3: "
			<< +pow(abs(complex<double>(W[f_index][0][1][0], -W[f_index][0][1][1])
			*complex<double>(H[f_index][1][0][0], H[f_index][1][0][1])), 2)*sigma_rk2
			+ pow(abs(complex<double>(W[f_index][0][1][0], -W[f_index][0][1][1])
			*complex<double>(H[f_index][1][1][0], H[f_index][1][1][1])), 2)*sigma_rk2
			+ pow(abs(complex<double>(W[f_index][0][1][0], -W[f_index][0][1][1])
			*complex<double>(H[f_index][1][2][0], H[f_index][1][2][1])), 2)*sigma_rk2 << endl;
			cout << "W2*: "
			<< complex<double>(W[f_index][0][0][0], -W[f_index][0][0][1]) << endl;
			cout << "W3*: "
			<< complex<double>(W[f_index][0][1][0], -W[f_index][0][1][1]) << endl;
			cout << "Sigma_nu2: " << sigma_nu2 << endl;
			cout << "Sigma_nu: " << sigma_nu << endl;
			cout << "Sigma_rk2: " << sigma_rk2 << endl;
			cout << "CNR: " << CNR << endl;
			cout << "Sigma_m121: " << sigma_m121 << endl;
			cout << "Sigma_m122: " << sigma_m122 << endl;
			cout << "a_1: " << a_1 << endl;
			cout << "a_2: " << a_2 << endl;
			cout << "a_3: " << a_3 << endl;
			cout << "a_4: " << a_4 << endl;
			cout << "a_1 + a_3: " << a_1 + a_3 << endl << "a_2 + a_4: " << a_2 + a_4 << endl;
			getchar();																		   */
		theo_BER = 0; //x_max = 0; x_min = 0;
		for (sb2hat = 0; sb2hat < pow(Nb, 2); sb2hat++)
		{
			for (sb2 = 0; sb2 < pow(Nb, 2); sb2++)
			{
				for (sb1 = 0; sb1 < pow(Nb, 2); sb1++)
				{
					for (sb1hat = 0; sb1hat < pow(Nb, 2); sb1hat++)
					{
						if (Symbols(sb1).real() != Symbols(sb1hat).real() || Symbols(sb1).imag() != Symbols(sb1hat).imag())
						{
							DeltaSb1 = (Symbols(sb1) - Symbols(sb1hat));
							DeltaSb2 = (Symbols(sb2) - Symbols(sb2hat));
							/* calculation of denominator */
							denominator = 2 * abs(a_1*DeltaSb1 + a_2*DeltaSb2)*sigma_nu;
							/* calculation of x */
							/*x = 2*(((a_1 + a_3)*Symbols(sb1) + (a_2 + a_4)*Symbols(sb2))*conj(-a_1*DeltaSb1 - a_2*DeltaSb2)).real();
							x += pow(abs(a_1*DeltaSb1 + a_2*DeltaSb2), 2);*/
							x = pow(abs(a_1*Symbols(sb1hat) + a_2*Symbols(sb2hat) + a_3*Symbols(sb1) + a_4*Symbols(sb2)), 2)
								- pow(abs((a_1 + a_3)*Symbols(sb1) + (a_2 + a_4)*Symbols(sb2)), 2);
							//x = abs(x);
							/* calculation of BER */
							Ne = (Symbols(sb1).real() != Symbols(sb1hat).real() && Symbols(sb1).imag() != Symbols(sb1hat).imag() ? 2 : 1);
							Pe = 0.5*erfc(x / denominator);
							theo_BER += Ne*Pe;

							/*cout << "Sb1: " << Symbols(sb1) << endl;
							cout << "Sb2: " << Symbols(sb2) << endl;
							cout << "Sb1hat: " << Symbols(sb1hat) << endl;
							cout << "Sb2hat: " << Symbols(sb2hat) << endl;
							cout << "Deltasb1: " << (Symbols(sb1) - Symbols(sb1hat)) << endl;
							cout << "Deltasb2: " << (Symbols(sb2) - Symbols(sb2hat)) << endl;
							cout << "denominator: " << denominator << endl;
							cout << "x = " << pow(abs(a_1*Symbols(sb1hat) + a_2*Symbols(sb2hat) + a_3*Symbols(sb1) + a_4*Symbols(sb2)), 2) << " - " << pow(abs((a_1 + a_3)*Symbols(sb1) + (a_2 + a_4)*Symbols(sb2)), 2) << endl;
							cout << "x/denom: " << x << " / " << denominator << endl;
							cout << "erfc( " << x / denominator << " ): " << erfc(x / denominator) << endl;
							cout << "number of erroneous bits: " << (abs(Symbols(sb1) - Symbols(sb1hat)) == 0
							? 0 : (abs(Symbols(sb1) - Symbols(sb1hat)) == 2 * OneBySqrt2 ? 1 : 2)) << endl;
							cout << "theo_ber: " << theo_BER << endl;
							cout << "residueal sb3: "
							<< pow(abs(complex<double>(W[f_index][0][0][0], -W[f_index][0][0][1])
							*complex<double>(H[f_index][0][2][0], H[f_index][0][2][1])
							+ complex<double>(W[f_index][0][1][0], -W[f_index][0][1][1])
							*complex<double>(H[f_index][1][2][0], H[f_index][1][2][1])), 2) << endl;
							getchar();	*/
							//if (x > x_max)x_max = x; if (x < x_min)x_min = x;
						}
					}
				}
			}
		}
		/* calculation of average BER over N*/
		theo_BER /= Nb;
		theo_BER_N = ((double)f_index / ((double)f_index + 1))*theo_BER_N + theo_BER / ((double)f_index + 1);
		//cout << "theoretical ber: " << theoretical_BER << endl;
		//cout << "x_max, x_min: " << x_max << ", " << x_min << endl;
		//getchar();
	}
	theoretical_BER = ((double)loop / ((double)loop + 1))*theoretical_BER + theo_BER_N / ((double)loop + 1);
	//getchar();
	cout << "Eb/N0=" << (int)CNR - 3 * MODULATION / 2 << "::loop=" << loop << ":::Theoretical BER =" << theoretical_BER << endl;

	/*if ((err = fopen_s(&theoretical_ber_data, "theoretical_ber_data.dat", "a")) != 0) { printf("The file outputlog file was not opened\n");	getchar();	exit(1); };

	file_time = currentDateTime();
	fprintf(theoretical_ber_data, "# %s ==  ", file_time.c_str());

	printf("Eb/N0=%2d::loop=%5d::Theoretical BER = %1.2e\n", (int)CNR - 3 * MODULATION / 2, loop, theoretical_BER);
	fprintf(theoretical_ber_data, "Eb/N0=%2d::loop=%5d::Theoretical BER = %1.2e\n", (int)CNR - 3 * MODULATION / 2, loop, theoretical_BER);
	fclose(theoretical_ber_data);*/
	if (loop == maxloop - 1) {
		cout << "Eb/N0=" << (int)CNR - 3 * MODULATION / 2 << ":::Average Theoretical BER =" << theoretical_BER << endl;

		if ((err = fopen_s(&theoretical_ber_data, "theoretical_BER.dat", "a")) != 0) { printf("The file Theoretical_BER.dat file was not opened\n");	getchar();	exit(1); };

		fprintf(theoretical_ber_data, "%1d\t%e\n", (int)CNR - 3 * MODULATION / 2, theoretical_BER);
		fclose(theoretical_ber_data);
	}
}
