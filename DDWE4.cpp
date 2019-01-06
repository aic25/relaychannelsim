#include "const.h"		 

extern double EbNo, CNR;
extern int loop;
extern double channel_gain[N][2], p_in[N];
extern errno_t err;	/* Error vector */
extern double P_for_DDWE[N][CL][CL][2], j[CL][2];
extern double OutPflag;
double OutPflag2 = 1;

const std::string currentDateTime() {
	// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	// Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
	// for more information about date/time format
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

	return buf;
}

void DDWE_Update(
	int(*Decision),
	int R_index,
	double(*received_vector)[2],
	double(*received_vector_wonoise)[2],
	double(*reference_symbol),
	double(*old_weights)[2],
	double(*new_weights)[2],
	double(*H)[BS][2],
	double(*e)[NUMBER_OF_STEPS][2],
	double(*e_w)[NUMBER_OF_STEPS][2],
	double(*e_w3)[NUMBER_OF_STEPS][2],
	double(*reference_symbol_bs3),
	int f_index,
	int index){

#ifndef DDWE_VARIABLES
	static double
		pf_findex,
		pf_counter;
	Matrix3cf TempoH;
	Vector3f EigValue_Differences;
	ComplexEigenSolver<Matrix3cf> TempoHEig;
	double
		DeltaW,
		power[2],
		gH[CL][2],
		g_i[BS][2],
		power_i[2],
		g_iH[BS][2],
		no = pow(10.0, -CNR / 10.0) / CODING_RATE,
		bir[2] = { 1, 0 },
		sifir[CL][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } },
		ff_inverse[2] = { FF_INVERSE, 0 },
		ff[2] = { FORGETTING_FACTOR, 0 },
		temp[2],
		temp2[2],
		temp3[2],
		received_symbol[2],
		R[BS][2],
		tempvec[CL][2],
		tempvech[CL][2],
		tempvec2[CL][2],
		tempvec3[CL][2],
		oldw_store[CL][2],
		x[CL][2],
		xh[CL][2],
		lpx[CL][2],
		lP_for_DDWE[N][CL][CL][2],
		tempmat[CL][CL][2],
		tempmat2[CL][CL][2],
		s_qpsk[4][2] = { { OneBySqrt2, OneBySqrt2 }, { OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, -OneBySqrt2 }, { -OneBySqrt2, OneBySqrt2 } };
	int s3_index, s1_index, i;
#endif // !DDWE_VARIABLES						   

#ifndef UPDATE_WEIGHTS
	for ( i = 0; i < TS; i++)
	{
		x[i][0] = received_vector[i][0];
		x[i][1] = received_vector[i][1];
	}
	x[CL - 1][0] = reference_symbol_bs3[0];
	x[CL - 1][1] = reference_symbol_bs3[1];


	VectorConjugate_CL(x, xh);
	ScalarMulMatrixToMatrix_CLxCL(ff_inverse, P_for_DDWE[f_index], lP_for_DDWE[f_index]);
	MatrixMulVectorToVector_CLxCL(lP_for_DDWE[f_index], x, lpx);	// nom
	VectorMulVectorToScalar_CL(xh, lpx, temp2);
	temp2[0] += bir[0];	// denom
	ScalarDivScalarToScalar(bir, temp2, temp3);
	ScalarMulVectorToVector_CL(temp3, lpx, j);	// eq. 1

	VectorEqualVector_CL(old_weights, oldw_store);

	VectorConjugate_CL(old_weights, tempvech);
	VectorMulVectorToScalar_CL(x, tempvech, received_symbol);
	//VectorMulVectorToScalar_CL(xh, old_weights, received_symbol);
	/*temp[0] = 0.0;
	temp[1] = 0.0;
	for (ts = 0; ts < CL; ts++)
	{
	temp[0] += W[f_index][0][ts][0] * tempvec[ts][0] - W[f_index][0][ts][1] * tempvec[ts][1];
	temp[1] += W[f_index][0][ts][0] * tempvec[ts][1] + W[f_index][0][ts][1] * tempvec[ts][0];
	} */

	e[f_index][index][0] = reference_symbol[0] - received_symbol[0];
	e[f_index][index][1] = -(reference_symbol[1] - received_symbol[1]);

#ifdef W_OUTLOG
	if ((sqr(e[f_index][index][0]) + sqr(e[f_index][index][1]) > W_OUTPUT_TRESHOLD) || ((pf_counter > 0) && (pf_findex == f_index))){
		if (pf_counter == 0) { pf_findex = f_index; pf_counter = 5; cout << "if1:" << pf_counter; }
		if (pf_counter > 0) { pf_counter--; cout << "if2" << pf_counter; }
		cout << endl;
		cout << "Counter: " << pf_counter << endl;
		cout << "f_index: " << f_index << "\tindex: " << index << endl;
		cout << received_symbol[0] << "\t" << received_symbol[1] << endl;
		cout << reference_symbol[0] << "\t" << reference_symbol[1] << endl;
		cout << reference_symbol_bs3[0] << "\t" << reference_symbol_bs3[1] << endl;
		cout << "Received vector:" << endl;
		cout << received_vector[0][0] << "\t" << received_vector[0][1] << "\tPower: " << sqr(received_vector[0][1]) + sqr(received_vector[0][0]) << endl;
		cout << received_vector[1][0] << "\t" << received_vector[1][1] << "\tPower: " << sqr(received_vector[1][1]) + sqr(received_vector[1][0]) << endl;
		cout << "Old weights:" << endl;
		cout << old_weights[0][0] << "\t" << old_weights[0][1] << endl;
		cout << old_weights[1][0] << "\t" << old_weights[1][1] << endl;
		cout << old_weights[2][0] << "\t" << old_weights[2][1] << endl;
		VectorMulVectorToScalarH_CL(old_weights, old_weights, temp3);
		cout << "Weight vector power: " << temp3[0] << endl;
	}
#endif // W_OUTLOG	  

	//if (sqr(e[f_index][index][0]) + sqr(e[f_index][index][1]) < 1)
	//{
	ScalarMulVectorToVector_CL(e[f_index][index], j, tempvec);
	VectorAddVectorToVector_CL(old_weights, tempvec, tempvec2);
	VectorAddVectorToVector_CL(tempvec2, sifir, new_weights);

	VectorMulVectorToScalarH_CL(new_weights, new_weights, e_w[f_index][index]);
	ScalarMulScalar_H(new_weights[2], new_weights[2], e_w3[f_index][index]);

	VectorMulMatrixToVector_CLxCL(xh, lP_for_DDWE[f_index], tempvec);
	VectorMulVectorToMatrix_CLxCL(j, tempvec, tempmat);
	MatrixSubMatrixToMatrix_CL(lP_for_DDWE[f_index], tempmat, P_for_DDWE[f_index]);
	//}

#ifdef RCV_VEC_AND_RCV_SIG_OUTPUT
	FILE *Rec_Vec;
	char  FName_RV[50] = "DDWEvalues_";
	char Ext[5] = ".dat";
	char cnr[10];
	_itoa(CNR, cnr, 10);
	std::strcat(FName_RV, cnr);

	/* it is possible add loop and f_index to name and make the file for each of those */

	std::strcat(FName_RV, Ext);
	if ((err = fopen_s(&Rec_Vec, FName_RV, "a+")) != 0){
		printf("The file Rec_Vec.dat file was not opened\n");
		getchar();
		exit(1);
	};					   	
	
	
	






























	
	/*if (OutPflag2 == 1){
		fprintf(Rec_Vec, "
		#Step\t
		#e(det-dec)1\t
		#e(det-dec)2\t
		#e(det-dec)1+2\t
		#IW_n-W_oI\t
		#x[0]\t
		#x[1]\t
		#IxI\t
		#d[0]\t
		#d[1]\t
		#IdI\t
		#IeI\t
		#Iw1I\t
		#Iw2I\t
		#Iw3I\t
		#IwI\t	 
		#Iw1_oI\t
		#Iw2_oI\t
		#Iw3_oI\t
		#Iw_oI\t
		#eigval1\t
		#eigval2\t
		#eigval3\t
		#sumeigval\t
		#eigmaxdiff\t
		#|H_1^(1)|\t	  
		#|H_2^(1)|\t
		#|H_3^(1)|\t
		#|H_1^(2)|\t
		#|H_2^(2)|\t
		#|H_3^(2)|\t
		#x[0][0]\t
		#x[0][1]\t
		#x[1][0]\t
		#x[1][1]\t	   
		#xwonoise[0][0]\t
		#xwonoise[0][1]\t
		#xwonoise[1][0]\t
		#xwonoise[1][1]\t
		\n");
		OutPflag2 = 0;
	} */
	if (sqrt(sqr(e[f_index][index][0]) + sqr(e[f_index][index][1])) > POUT_THOLD)
	{
		/*TempoH(0, 0).real(H[0][0][0]); TempoH(0, 0).imag(H[0][0][1]);
		TempoH(0, 1).real(H[0][1][0]); TempoH(0, 1).imag(H[0][1][1]);
		TempoH(0, 2).real(H[0][2][0]); TempoH(0, 2).imag(H[0][2][1]);
		TempoH(1, 0).real(H[1][0][0]); TempoH(1, 0).imag(H[1][0][1]);
		TempoH(1, 1).real(H[1][1][0]); TempoH(1, 1).imag(H[1][1][1]);
		TempoH(1, 2).real(H[1][2][0]); TempoH(1, 2).imag(H[1][2][1]);
		TempoH(2, 0).real(0); TempoH(2, 0).imag(0);
		TempoH(2, 1).real(1); TempoH(2, 1).imag(0);
		TempoH(2, 2).real(0); TempoH(2, 2).imag(0);
		TempoH = TempoH*TempoH.adjoint();
		TempoHEig.compute(TempoH);
		EigValue_Differences(0) = abs(TempoHEig.eigenvalues()[0].real() - TempoHEig.eigenvalues()[1].real());
		EigValue_Differences(1) = abs(TempoHEig.eigenvalues()[0].real() - TempoHEig.eigenvalues()[2].real());
		EigValue_Differences(2) = abs(TempoHEig.eigenvalues()[1].real() - TempoHEig.eigenvalues()[2].real());
		
		cout << "Tempo: " << endl << TempoH << endl;
		cout << "Tempo: " << endl << TempoH << endl;
		cout << "Tempo Eigenvals: " << endl << TempoHEig.eigenvalues() << endl;
		//cout << "Channel matrix reobtained: " << endl << TempoHEig.eigenvectors()*TempoHEig.eigenvalues().asDiagonal()*TempoHEig.eigenvectors().inverse() << endl;
		cout << "eig trace: " << TempoHEig.eigenvalues().sum() << endl;
		cout << "Eigval diff: " << endl << EigValue_Differences << endl;
		cout << "Max diff: " << endl << EigValue_Differences.maxCoeff() << endl;
		getchar();	  																  
		*/	
		DeltaW = MagnitudeOfDifference_CL(oldw_store, new_weights);
		fprintf(Rec_Vec, "%d(%d)\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
			index,f_index,
			(sqrt(sqr(reference_symbol[0] - s_qpsk[Decision[0]][0]) + sqr(reference_symbol[1] - s_qpsk[Decision[0]][1]))),
			(sqrt(sqr(reference_symbol_bs3[0] - s_qpsk[Decision[1]][0]) + sqr(reference_symbol_bs3[1] - s_qpsk[Decision[1]][1]))),
			(sqrt(sqr(reference_symbol[0] - s_qpsk[Decision[0]][0]) + sqr(reference_symbol[1] - s_qpsk[Decision[0]][1]))) + (sqrt(sqr(reference_symbol_bs3[0] - s_qpsk[Decision[1]][0]) + sqr(reference_symbol_bs3[1] - s_qpsk[Decision[1]][1]))),
			DeltaW,
			sqrt(sqr(x[0][0]) + sqr(x[0][1])),
			sqrt(sqr(x[1][0]) + sqr(x[1][1])),
			sqrt(sqr(x[0][0]) + sqr(x[0][1]) + sqr(x[1][0]) + sqr(x[1][1]) + sqr(x[2][0]) + sqr(x[2][1])),
			received_symbol[0],
			received_symbol[1],
			sqrt(sqr(received_symbol[0]) + sqr(received_symbol[1])),
			sqrt(sqr(e[f_index][index][0]) + sqr(e[f_index][index][1])),
			sqrt(sqr(new_weights[0][0]) + sqr(new_weights[0][1])),
			sqrt(sqr(new_weights[1][0]) + sqr(new_weights[1][1])),
			sqrt(sqr(new_weights[2][0]) + sqr(new_weights[2][1])),
			sqrt(sqr(new_weights[0][0]) + sqr(new_weights[0][1]) + sqr(new_weights[1][0]) + sqr(new_weights[1][1]) + sqr(new_weights[2][0]) + sqr(new_weights[2][1])),
			sqrt(sqr(oldw_store[0][0]) + sqr(oldw_store[0][1])),
			sqrt(sqr(oldw_store[1][0]) + sqr(oldw_store[1][1])),
			sqrt(sqr(oldw_store[2][0]) + sqr(oldw_store[2][1])),
			sqrt(sqr(oldw_store[0][0]) + sqr(oldw_store[0][1]) + sqr(oldw_store[1][0]) + sqr(oldw_store[1][1]) + sqr(oldw_store[2][0]) + sqr(oldw_store[2][1])),
			//TempoHEig.eigenvalues()[0].real(),
			//TempoHEig.eigenvalues()[1].real(),
			//TempoHEig.eigenvalues()[2].real(),
			//TempoHEig.eigenvalues().sum().real(),
			//EigValue_Differences.maxCoeff(),		   
			0.0,
			0.0,
			0.0,
			0.0,
			0.0,
			sqrt(sqr(H[0][0][0]) + sqr(H[0][0][1])),
			sqrt(sqr(H[0][1][0]) + sqr(H[0][1][1])),
			sqrt(sqr(H[0][2][0]) + sqr(H[0][2][1])),
			sqrt(sqr(H[1][0][0]) + sqr(H[1][0][1])),
			sqrt(sqr(H[1][1][0]) + sqr(H[1][1][1])),
			sqrt(sqr(H[1][2][0]) + sqr(H[1][2][1])),	 
			x[0][0],
			x[0][1],
			x[1][0],
			x[1][1],
			received_vector_wonoise[0][0],
			received_vector_wonoise[0][1],
			received_vector_wonoise[1][0],
			received_vector_wonoise[1][1]
			);

	}	fclose(Rec_Vec);
#endif // RCV_VEC_AND_RCV_SIG_OUTPUT  

#ifdef W_OUTLOG
	if ((sqr(e[f_index][index][0]) + sqr(e[f_index][index][1]) > W_OUTPUT_TRESHOLD) || ((pf_counter > 0) && (pf_findex == f_index))){
		cout << "New weights:" << endl;
		cout << new_weights[0][0] << "\t" << new_weights[0][1] << endl;
		cout << new_weights[1][0] << "\t" << new_weights[1][1] << endl;
		cout << new_weights[2][0] << "\t" << new_weights[2][1] << endl;
		VectorMulVectorToScalarH_CL(new_weights, new_weights, temp3);
		cout << "Weight vector power: " << temp3[0] << endl;
		OutPflag = loop;
		getchar();
	}
#endif // W_OUTLOG

#endif // !UPDATE_WEIGHTS

#ifndef UPDATE_POWER
#ifdef SD
	VectorMulMatrixToVector_TSxTS(new_weights, H, R);
	g_i[0][0] = 0.0; g_i[0][1] = 0.0;
	for (size_t i = 1; i < M; i++){g_i[i][0] = R[i][0]; g_i[i][1] = R[i][1];}
	VectorConjugate_M(g_i, g_iH);
	VectorMulVectorToScalar_M(g_iH, g_i, power_i);	/* power_i = g_i^H x g_i */
	VectorConjugate_TS(new_weights, gH);
	VectorMulVectorToScalar_TS(gH, new_weights, power);
	channel_gain[f_index][0] = R[0][0];
	channel_gain[f_index][1] = -R[0][1];
	p_in[f_index] = no * power[0] + power_i[0];	/* Incoming signal power. */
#endif // SD						  
#endif // !POWER_UPDATE
}	/* DDWE_Update :: end */