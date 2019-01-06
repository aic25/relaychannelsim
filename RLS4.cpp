#include "const.h"
extern double CNR;
extern double channel_gain[N][2], p_in[N];
double EigSort[4];
int icerik[4];
double P_for_DDWE[N][CL][CL][2], j[CL][2];
Matrix4cf XCorr[N];
//float minEigVal;

extern errno_t err;	/* Error vector */

void RLS_initial_weigth(
	int f_index,
	int receiver_index,
	double(*DFT_signal)[TS][N][2],
	double(*reference_signal)[TS][N][2],
	double(*W)[4][4][2],
	double(*H)[TS][BS][2],
	double(*wh_sequence)[Np][N],
	double(*e)[NUMBER_OF_STEPS][2],
	double(*e_w)[NUMBER_OF_STEPS][2],
	double(*e_w3)[NUMBER_OF_STEPS][2]){

#ifndef RLS_VARIABLES
	int	i, j, np, v, ts, testi, temp;
	double n;
	Matrix4cf XXH;
	Matrix4cf Test;
	Vector4cf X;
	ComplexEigenSolver<Matrix4cf> Eigs;
	std::ptrdiff_t a, b;
	//std::complex<float> minEigVal;
	double
		sigma = 0.01,
		no,
		bir[2] = { 1, 0 },
		sifir[CL][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } },
		ff_inverse[2] = { FF_INVERSE, 0 },
		ff[2] = { FORGETTING_FACTOR, 0 },
		x[4][2],
		xh[4][2],
		xxH[4][4][2];
#endif // !RLS_VARIABLES

#ifndef RLS_INITIALIZATION
	for (i = 0; i < 4; i++)	{
		for (j = 0; j < 4; j++)	{
			XCorr[f_index](i, j).real(0);
			XCorr[f_index](i, j).imag(0);
		}
	}
#endif // !RLS_INITIALIZATION

	/* -------------------------------------------- */
	no = pow(10.0, -CNR / 10.0) / CODING_RATE;	/* Noise power. Eb = 1 */

	for (n = 1; n < PREAMBLE_LENGTH + 1; n++){
		np = n - 1;
		/*X[0].real(DFT_signal[np][0][f_index][0]);	X[0].imag(DFT_signal[np][0][f_index][1]);
		X[1].real(DFT_signal[np][1][f_index][0]);	X[1].imag(DFT_signal[np][1][f_index][1]);
		X[2].real(wh_sequence[(receiver_index + 1) % 2][np][f_index]);	X[2].imag(0);
		X[3].real(wh_sequence[receiver_index][np][f_index]);	X[3].imag(0);
		cout << "X(in): " << endl << X << endl;		*/
		X[0] = complex<double>(DFT_signal[np][0][f_index][0], DFT_signal[np][0][f_index][1]);
		X[1] = complex<double>(DFT_signal[np][1][f_index][0], DFT_signal[np][1][f_index][1]);
		X[2] = complex<double>(wh_sequence[(receiver_index + 1) % 2][np][f_index], 0);
		X[3] = complex<double>(wh_sequence[receiver_index][np][f_index], 0);
		//cout << "X=(out): " << endl << X << endl; getchar();

		XXH = X*X.adjoint();
		/*for (i = 0; i < 4; i++)
			for (j = 0; j < 4; j++){
			XXH(i, j).real(X(i).real()*X(j).real() + X(i).imag()*X(j).imag());
			XXH(i, j).imag(X(i).imag()*X(j).real() - X(i).real()*X(j).imag());
			}

			if (f_index == 1) {
			cout << "X:" << endl << X << endl
			//<< "X^H: " << endl << X.adjoint() << endl
			//<< "1/5:" << 1 / 5 << " 1/test(4): " << 1 / test << " 1/test2(3): " << 1 / test2 << " 1/testi(6) " << 1 / testi << endl
			<< "F: " << f_index << " " << "step: " << n << " 1/step= " << 1 / n << (1 / n) << endl
			//<< " XCorr[f]: " << endl << XCorr[f_index] << endl
			<< "XXH: " << endl << XXH << endl
			<< "1/nXXH:" << (1 / n) << endl << (1 / n)*XXH << endl
			//<< "(1/n)XX.H: " << endl << (1 / n)*X*X.adjoint() << endl
			//<< "Test: " << endl << Test << endl
			//<< "?: " << X[0] * X[0] << endl
			//<< "n-1/nXCorr: " << endl << ((n - 1) / n)*XCorr[f_index] << endl
			//<< "and: "
			<< endl;
			getchar();
			}	  */
		XCorr[f_index] = (1 / n)*XXH + ((n - 1) / n)*XCorr[f_index];
		//if (f_index == 1) { cout << XCorr[f_index] << endl; getchar(); }
		//if (f_index == 1){	cout << XXH << endl << XCorr << endl; }	  
		//if (f_index == 5){ cout << "XCorr: " << endl << XCorr[f_index] << endl; getchar(); }
	}	/* end :: n */
	/*XCorr[f_index](2, 2).real(1); XCorr[f_index](2, 2).imag(0);
	XCorr[f_index](2, 3).real(0); XCorr[f_index](2, 2).imag(0);
	XCorr[f_index](3, 2).real(0); XCorr[f_index](2, 2).imag(0);
	XCorr[f_index](3, 3).real(1); XCorr[f_index](2, 2).imag(0);	   */
	/*XCorr[f_index](2, 2) = complex<double>(1, 0);
	XCorr[f_index](2, 3) = complex<double>(0, 0);
	XCorr[f_index](3, 2) = complex<double>(0, 0);
	XCorr[f_index](3, 3) = complex<double>(1, 0);	  */
	//if (f_index == 5){cout << "XCorr: " << endl << XCorr[f_index] << endl; getchar();	}
	Eigs.compute(XCorr[f_index]);
	//minEigVal = Eigs.eigenvalues().real().minCoeff(&a);
	EigSort[0] = Eigs.eigenvalues().real()[0];
	EigSort[1] = Eigs.eigenvalues().real()[1];
	EigSort[2] = Eigs.eigenvalues().real()[2];
	EigSort[3] = Eigs.eigenvalues().real()[3];
	sirala(EigSort, 4, icerik);
	EigSort[0] = Eigs.eigenvalues().real()[0];
	EigSort[1] = Eigs.eigenvalues().real()[1];
	EigSort[2] = Eigs.eigenvalues().real()[2];
	EigSort[3] = Eigs.eigenvalues().real()[3];

	/*cout << Eigs.eigenvalues().real() << endl
		<< EigSort[0] << endl
		<< EigSort[1] << endl
		<< EigSort[2] << endl
		<< EigSort[3] << endl
		<< icerik[0] << endl
		<< icerik[1] << endl
		<< icerik[2] << endl
		<< icerik[3] << endl;
		getchar();*/

	/*if (f_index == 1) {
		cout << minEigVal << " " << a << endl
		<< "XCorr: " << endl
		<< XCorr[f_index] << endl
		<< "Eigvals: " << endl
		<< Eigs.eigenvalues() << endl
		<< "Eigvecs: " << endl
		<< Eigs.eigenvectors() << endl
		<< "Ax1= " << endl
		<< XCorr[f_index]*Eigs.eigenvectors().col(0) << endl
		<< "lamda1x1= " << endl
		<< Eigs.eigenvalues()[0] * Eigs.eigenvectors().col(0) << endl
		<< "Ax2= " << endl
		<< XCorr[f_index] * Eigs.eigenvectors().col(1) << endl
		<< "lamda2x2= " << endl
		<< Eigs.eigenvalues()[1] * Eigs.eigenvectors().col(1) << endl
		<< "Ax3= " << endl
		<< XCorr[f_index] * Eigs.eigenvectors().col(2) << endl
		<< "lamda3x3= " << endl
		<< Eigs.eigenvalues()[2] * Eigs.eigenvectors().col(2) << endl
		<< "Ax4= " << endl
		<< XCorr[f_index] * Eigs.eigenvectors().col(3) << endl
		<< "lamda4x4= " << endl
		<< Eigs.eigenvalues()[3] * Eigs.eigenvectors().col(3) << endl
		;
		getchar(); }	*/
	/*cout << Eigs.eigenvectors().col(a) << endl
		<< Eigs.eigenvectors().col(a)[0] << endl
		<< Eigs.eigenvectors().col(a)[0].real() << endl; */

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

	/*cout << "Normalized eig vec: " << endl
		<< Eigs.eigenvectors().col(0).normalized() << endl
		<< "EV_N^H*EV_N: " << endl
		<< Eigs.eigenvectors().col(0).adjoint()*Eigs.eigenvectors().col(0) << endl
		<< "Not normalized: " << endl
		<< Eigs.eigenvectors().col(0) << endl
		<< "EV^H*EV: " << endl
		<< Eigs.eigenvectors().col(0).adjoint()*Eigs.eigenvectors().col(0) << endl;
		getchar();	 */
	/*cout << "Eigenvalues:" << endl
		<< Eigs.eigenvalues() << endl
		<< "Test eigval output:" << endl
		<< minEigVal << endl
		<< "Corresponding EigVec: " << endl
		<< Eigs.eigenvectors().col(a) << endl
		<< "Transferred C matrix output: " << endl
		<< W[f_index][0][0][0] << " " << W[f_index][0][0][1] << endl
		<< W[f_index][0][1][0] << " " << W[f_index][0][1][1] << endl
		<< W[f_index][0][2][0] << " " << W[f_index][0][2][1] << endl
		<< W[f_index][0][3][0] << " " << W[f_index][0][3][1] << endl;

		getchar();															 */
}
