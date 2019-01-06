#include "const-function-declerations.h"  
extern double channel_gain[N][2], p_in[N], wh_sequence[BS][Np][N][2], CNR;
double EigSort[4];
int icerik[4];
Matrix4cf XCorr[N];
extern errno_t err;	/* Error vector */

void RLS_initial_weigth(
	int f_index,
	int receiver_index,
	double(*DFT_signal)[TS][N][2],
	double(*reference_signal)[TS][N][2],
	double(*W)[4][4][2],
	double(*H)[TS][BS][2]){
											  
	double W_PL4[4][2], W_PL64[4][2];
	double temp[2];

	int	i, j, np, v, ts;
	double n, eigval_store[2];
	Matrix4cf XXH;
	Matrix4cf Temp4cfM;
	Matrix4cf Test;
	Vector4cf X;
	ComplexEigenSolver<Matrix4cf> Eigs;
	std::ptrdiff_t a, b;
	double
		bir[2] = { 1, 0 },
		sifir[CL][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } },
		ff_inverse[2] = { FF_INVERSE, 0 },
		ff[2] = { FORGETTING_FACTOR, 0 };

	for (i = 0; i < 4; i++)	{
		for (j = 0; j < 4; j++)	{
			XCorr[f_index](i, j).real(0);
			XCorr[f_index](i, j).imag(0);
		}
	}
	/* -------------------------------------------- */
	for (n = 1; n < PREAMBLE_LENGTH + 1; n++){
		np = n - 1;
		X[0] = complex<double>(DFT_signal[np][0][f_index][0], DFT_signal[np][0][f_index][1]);
		X[1] = complex<double>(DFT_signal[np][1][f_index][0], DFT_signal[np][1][f_index][1]);
		X[2] = complex<double>(wh_sequence[receiver_index][np][f_index][0], wh_sequence[receiver_index][np][f_index][1]);
		X[3] = complex<double>(wh_sequence[(receiver_index + 1) % 2][np][f_index][0], wh_sequence[(receiver_index + 1) % 2][np][f_index][1]);
		XXH = X*X.adjoint();
		Temp4cfM = XCorr[f_index];
		XCorr[f_index] = (1 / n)*XXH + ((n - 1) / n)*Temp4cfM;
		/*cout << "n: " << n << endl;
		cout << "xxh: " << endl << XXH << endl;
		cout << "1/n xxh: " << endl << (1 / n)*XXH << endl;
		cout << "xcorr: " << endl << XCorr[f_index] << endl;	getchar();		*/
		//if (n  4){
			Eigs.compute(XCorr[f_index]);
			/*EigSort[0] = Eigs.eigenvalues().real()[0];
			EigSort[1] = Eigs.eigenvalues().real()[1];
			EigSort[2] = Eigs.eigenvalues().real()[2];
			EigSort[3] = Eigs.eigenvalues().real()[3];
			sirala(EigSort, 4, icerik);
			EigSort[0] = Eigs.eigenvalues().real()[0];
			EigSort[1] = Eigs.eigenvalues().real()[1];
			EigSort[2] = Eigs.eigenvalues().real()[2];
			EigSort[3] = Eigs.eigenvalues().real()[3];
			W_PL4[0][0] = Eigs.eigenvectors().col(icerik[0])[0].real();
			W_PL4[0][1] = Eigs.eigenvectors().col(icerik[0])[0].imag();
			W_PL4[1][0] = Eigs.eigenvectors().col(icerik[0])[1].real();
			W_PL4[1][1] = Eigs.eigenvectors().col(icerik[0])[1].imag();
			W_PL4[2][0] = Eigs.eigenvectors().col(icerik[0])[2].real();
			W_PL4[2][1] = Eigs.eigenvectors().col(icerik[0])[2].imag();
			W_PL4[3][0] = Eigs.eigenvectors().col(icerik[0])[3].real();
			W_PL4[3][1] = Eigs.eigenvectors().col(icerik[0])[3].imag();
		}*/				  
	}	/* end :: n */
	//do {
		Eigs.compute(XCorr[f_index]);
		EigSort[0] = Eigs.eigenvalues().real()[0];
		EigSort[1] = Eigs.eigenvalues().real()[1];
		EigSort[2] = Eigs.eigenvalues().real()[2];
		EigSort[3] = Eigs.eigenvalues().real()[3];
		sirala(EigSort, 4, icerik);
		EigSort[0] = Eigs.eigenvalues().real()[0];
		EigSort[1] = Eigs.eigenvalues().real()[1];
		EigSort[2] = Eigs.eigenvalues().real()[2];
		EigSort[3] = Eigs.eigenvalues().real()[3];
		/*eigval_store[0] = EigSort[icerik[0]];
		Eigs.compute(XCorr[f_index]);
		EigSort[0] = Eigs.eigenvalues().real()[0];
		EigSort[1] = Eigs.eigenvalues().real()[1];
		EigSort[2] = Eigs.eigenvalues().real()[2];
		EigSort[3] = Eigs.eigenvalues().real()[3];
		sirala(EigSort, 4, icerik);
		EigSort[0] = Eigs.eigenvalues().real()[0];
		EigSort[1] = Eigs.eigenvalues().real()[1];
		EigSort[2] = Eigs.eigenvalues().real()[2];
		EigSort[3] = Eigs.eigenvalues().real()[3];
		eigval_store[1] = EigSort[icerik[0]];			*/
	//} while (eigval_store[0] < eigval_store[1]);


	/*W_PL64[0][0] = Eigs.eigenvectors().col(icerik[0])[0].real();
	W_PL64[0][1] = Eigs.eigenvectors().col(icerik[0])[0].imag();
	W_PL64[1][0] = Eigs.eigenvectors().col(icerik[0])[1].real();
	W_PL64[1][1] = Eigs.eigenvectors().col(icerik[0])[1].imag();
	W_PL64[2][0] = Eigs.eigenvectors().col(icerik[0])[2].real();
	W_PL64[2][1] = Eigs.eigenvectors().col(icerik[0])[2].imag();
	W_PL64[3][0] = Eigs.eigenvectors().col(icerik[0])[3].real();
	W_PL64[3][1] = Eigs.eigenvectors().col(icerik[0])[3].imag();

	temp[0] = sqr(W_PL4[0][0] - W_PL64[0][0]) +
		sqr(W_PL4[1][0] - W_PL64[1][0]) +
		sqr(W_PL4[2][0] - W_PL64[2][0]) +
		sqr(W_PL4[3][0] - W_PL64[3][0]);
	temp[1] = sqr(W_PL4[0][1] - W_PL64[0][1]) +
		sqr(W_PL4[1][1] - W_PL64[1][1]) +
		sqr(W_PL4[2][1] - W_PL64[2][1]) +
		sqr(W_PL4[3][1] - W_PL64[3][1]);
	if (sqrt(temp[0] + temp[1]) > 1){ 
		cout << "mag of difference: " << sqrt(temp[0] + temp[1]) << endl;
		cout << "w pl4: " << endl
			<< W_PL4[0][0] << "+j" << W_PL4[0][1] << " " << sqrt(sqr(W_PL4[0][0]) + sqr(W_PL4[0][1])) << endl
			<< W_PL4[1][0] << "+j" << W_PL4[1][1] << " " << sqrt(sqr(W_PL4[1][0]) + sqr(W_PL4[1][1])) << endl
			<< W_PL4[2][0] << "+j" << W_PL4[2][1] << " " << sqrt(sqr(W_PL4[2][0]) + sqr(W_PL4[2][1])) << endl
			<< W_PL4[3][0] << "+j" << W_PL4[3][1] << " " << sqrt(sqr(W_PL4[3][0]) + sqr(W_PL4[3][1])) << endl
			<< "w pl64: " << endl
			<< W_PL64[0][0] << "+j" << W_PL64[0][1] << " " << sqrt(sqr(W_PL64[0][0]) + sqr(W_PL64[0][1])) << endl
			<< W_PL64[1][0] << "+j" << W_PL64[1][1] << " " << sqrt(sqr(W_PL64[1][0]) + sqr(W_PL64[1][1])) << endl
			<< W_PL64[2][0] << "+j" << W_PL64[2][1] << " " << sqrt(sqr(W_PL64[2][0]) + sqr(W_PL64[2][1])) << endl
			<< W_PL64[3][0] << "+j" << W_PL64[3][1] << " " << sqrt(sqr(W_PL64[3][0]) + sqr(W_PL64[3][1])) << endl;
		getchar();
	}				   */

	W[f_index][receiver_index][0][0] = Eigs.eigenvectors().col(icerik[0])[0].real();
	W[f_index][receiver_index][0][1] = Eigs.eigenvectors().col(icerik[0])[0].imag();
	W[f_index][receiver_index][1][0] = Eigs.eigenvectors().col(icerik[0])[1].real();
	W[f_index][receiver_index][1][1] = Eigs.eigenvectors().col(icerik[0])[1].imag();
	W[f_index][receiver_index][2][0] = Eigs.eigenvectors().col(icerik[0])[2].real();
	W[f_index][receiver_index][2][1] = Eigs.eigenvectors().col(icerik[0])[2].imag();
	W[f_index][receiver_index][3][0] = Eigs.eigenvectors().col(icerik[0])[3].real();
	W[f_index][receiver_index][3][1] = Eigs.eigenvectors().col(icerik[0])[3].imag();

	W[f_index][(receiver_index + 1) % 2][0][0] = Eigs.eigenvectors().col(icerik[1])[0].real();
	W[f_index][(receiver_index + 1) % 2][0][1] = Eigs.eigenvectors().col(icerik[1])[0].imag();
	W[f_index][(receiver_index + 1) % 2][1][0] = Eigs.eigenvectors().col(icerik[1])[1].real();
	W[f_index][(receiver_index + 1) % 2][1][1] = Eigs.eigenvectors().col(icerik[1])[1].imag();
	W[f_index][(receiver_index + 1) % 2][2][0] = Eigs.eigenvectors().col(icerik[1])[2].real();
	W[f_index][(receiver_index + 1) % 2][2][1] = Eigs.eigenvectors().col(icerik[1])[2].imag();
	W[f_index][(receiver_index + 1) % 2][3][0] = Eigs.eigenvectors().col(icerik[1])[3].real();
	W[f_index][(receiver_index + 1) % 2][3][1] = Eigs.eigenvectors().col(icerik[1])[3].imag();

	W[f_index][2][0][0] = Eigs.eigenvectors().col(icerik[2])[0].real();
	W[f_index][2][0][1] = Eigs.eigenvectors().col(icerik[2])[0].imag();
	W[f_index][2][1][0] = Eigs.eigenvectors().col(icerik[2])[1].real();
	W[f_index][2][1][1] = Eigs.eigenvectors().col(icerik[2])[1].imag();
	W[f_index][2][2][0] = Eigs.eigenvectors().col(icerik[2])[2].real();
	W[f_index][2][2][1] = Eigs.eigenvectors().col(icerik[2])[2].imag();
	W[f_index][2][3][0] = Eigs.eigenvectors().col(icerik[2])[3].real();
	W[f_index][2][3][1] = Eigs.eigenvectors().col(icerik[2])[3].imag();

	W[f_index][3][0][0] = Eigs.eigenvectors().col(icerik[3])[0].real();
	W[f_index][3][0][1] = Eigs.eigenvectors().col(icerik[3])[0].imag();
	W[f_index][3][1][0] = Eigs.eigenvectors().col(icerik[3])[1].real();
	W[f_index][3][1][1] = Eigs.eigenvectors().col(icerik[3])[1].imag();
	W[f_index][3][2][0] = Eigs.eigenvectors().col(icerik[3])[2].real();
	W[f_index][3][2][1] = Eigs.eigenvectors().col(icerik[3])[2].imag();
	W[f_index][3][3][0] = Eigs.eigenvectors().col(icerik[3])[3].real();
	W[f_index][3][3][1] = Eigs.eigenvectors().col(icerik[3])[3].imag();

#ifdef SD
	MatrixMulMatrixToMatrix_BSxTS(W[f_index], H[f_index], R);
	VectorConjugate_TS(W[f_index][0], gH);
	VectorMulVectorToScalar_TS(gH, W[f_index][0], power);
	power_i[0] = 0;
	for (size_t i = 1; i < M; i++)
	{
		power_i[0] += sqr(R[0][i][0]) + sqr(R[0][i][1]);	/* 0.5=no */
	}
	channel_gain[f_index][0] = R[0][0][0];
	channel_gain[f_index][1] = -R[0][0][1];
	p_in[f_index] = (no)*power[0] + power_i[0];
#endif // SD	  
	/* end :: f_index */
}
