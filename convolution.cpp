#include "const.h"

void zerotap(const double Signal[/* SignalLen */][2], double SignalExtended[/* SignalLen * SignalDelay */][2], size_t SignalLen, size_t SignalDelay){
	size_t m,n;
	for (m = 0; m < SignalLen; m++){
		SignalExtended[m*SignalDelay][0] = Signal[m][0];
		SignalExtended[m*SignalDelay][1] = Signal[m][1];
		for (n = 1; n < SignalDelay; n++){
			SignalExtended[m*SignalDelay + n][0] = 0.0;
			SignalExtended[m*SignalDelay + n][1] = 0.0;
		}
	}
}
void convolve(const double Signal[/* SignalLen */][2], size_t SignalLen,
	const double Kernel[/* KernelLen */][2], size_t KernelLen,
	double Result[/* SignalLen + KernelLen - 1 */][2])
{
	size_t m, n, o;

	for (n = 0; n < SignalLen + KernelLen - 1; n++)
	{
		size_t kmin, kmax, k;

		Result[n][0] = 0.0;
		Result[n][1] = 0.0;

		kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0.0;
		kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

		for (k = kmin; k <= kmax; k++)
		{
			Result[n][0] += Signal[k][0] * Kernel[n - k][0] - Signal[k][1] * Kernel[n - k][1];
			Result[n][1] += Signal[k][1] * Kernel[n - k][0] + Signal[k][0] * Kernel[n - k][1];
		}
	}
}

void printSignal(const char* Name, double Signal[/* SignalLen */][2], size_t SignalLen)
{
	size_t i;

	for (i = 0; i < SignalLen; i++)
	{
		//printf("%s[%zu].real = %f,\t%s[%zu].imag = %f,\n", Name, i, Signal[i][0], Name, i, Signal[i][1]);
		cout << Name << "[" << i << "].real" << Signal[i][0] << ", "
			<< Name << "[" << i << "].imag" << Signal[i][1] << endl;
	}
	printf("\n");
	//getchar();
}

/*int main(void)
{
double signal[] = { 1, 1, 1, 1, 1 };
double kernel[] = { 1, 1, 1, 1, 1 };
double result[ELEMENT_COUNT(signal) + ELEMENT_COUNT(kernel) - 1];

convolve(signal, ELEMENT_COUNT(signal),
kernel, ELEMENT_COUNT(kernel),
result);

printSignal("signal", signal, ELEMENT_COUNT(signal));
printSignal("kernel", kernel, ELEMENT_COUNT(kernel));
printSignal("result", result, ELEMENT_COUNT(result));

return 0;
}*/
