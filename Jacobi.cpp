#include "const.h"					


double Jacobi_Complex(double(***a), double EPS, double(***v), double(**d)){
	// a is the main square matrix with size M_Size x M_Size. v for eigen vectors and d for eigen values
	// still only for a real symmetric matrix


	int i, j, ip, iq, nrot;
	double tresh, theta, tau, t, sm, s, h, g, c;
	double b[M_Size][2], z[M_Size][2];

	for (i = 0; i < M_Size; i++){
		d[i][0] = 0; d[i][1] = 0;
		b[i][0] = 0; b[i][1] = 0;
		z[i][0] = 0; z[i][1] = 0;
		for (j = 0; j < M_Size; j++){
			v[i][j][0] = 0;
			v[i][j][1] = 0;

		}
	}

	for (ip = 0; ip < M_Size; ip++){
		for (iq = 0; iq < M_Size; iq++)v[ip][iq] = 0;
		v[ip][ip][0] = 1.0;
	}
	for (ip = 0; ip < M_Size; ip++){
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	for (i = 1; i < 50; ip++){
		sm = 0.0;
		for (ip = 0; ip < M_Size - 1; ip++){
			for (iq = ip + 1; iq < M_Size; iq++)
				sm += sqrt(sqr(a[ip][ip][0]) + sqr(a[ip][ip][1]));
		}
		if (sm == 0.0){
			eigsrt(d, v);
			return 0;
		}
		else
			tresh = 0.0;
		for (ip = 0; ip < M_Size - 1; ip++){
			for (iq = ip + 1; iq < M_Size; iq++){
				g = 100.0*abs(a[ip][iq]);
				if (i>4 && g <= EPS*abs(d[ip]) && g <= EPS*abs(d[iq]))
					a[ip][iq] = 0.0;
				else if (abs(a[ip][iq])>tresh){
					h = d[iq] - d[ip];
					if (g <= EPS*abs(h))
						t = (a[ip][iq]) / h;
					else {
						theta = 0.5*h / (a[ip][iq]);
						t = 1.0 / (abs(theta) + sqrt(1.0 + sqr(theta)));
						if (theta < 0.0)t = -t;
					}
					c = 1.0 / sqrt(1 + sqr(t));
					s = t*c;
					tau = s / (1.0 + c);
					h = t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
					for (j = 0; j < ip; j++)
						rot(a, s, tau, j, ip, j, iq);
					for (j = ip + 1; j < iq; j++)
						rot(a, s, tau, ip, j, j, iq);
					for (j = iq + i; j < M_Size; j++)
						rot(a, s, tau, ip, j, iq, j);
					for (j = 0; j < M_Size; j++)
						rot(v, s, tau, j, ip, j, iq);
					++nrot;
				}
			}
		}
		for (ip = 0; ip < M_Size; ip++){
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	return nrot;
}

void Rot_Complex(double(**a), const double s, const double tau, const int i, const int j, const int k, const int l){
	double g = a[i][j];
	double h = a[k][l];
	a[i][j] = g - s*(h + g*tau);
	a[k][l] = h + s*(g - h*tau);
}
void Eigsrt_Complex(double(*d), double(**v)){
	int k;
	int n = M_Size;
	for (int i = 0; i < n - 1; i++){
		double p = d[k = i];
		for (int j = i; j < n; j++)
		if (d[j] >= p)p = d[k = j];
		if (k != i){
			d[k] = d[i];
			d[i] = p;
			if (v != NULL)
			for (int j = 0; j < n; j++){
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}

double Jacobi_Real(double(**a), double EPS, double(**v), double(*d)){
	// a is the main square matrix with size M_Size x M_Size. v for eigen vectors and d for eigen values
	// still only for a real symmetric matrix


	int i, j, ip, iq, nrot;
	double tresh, theta, tau, t, sm, s, h, g, c;
	double b[M_Size], z[M_Size];

	for (ip = 0; ip < M_Size; ip++){
		for (iq = 0; iq < M_Size; iq++)v[ip][iq] = 0;
		v[ip][ip] = 1.0;
	}
	for (ip = 0; ip < M_Size; ip++){
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	for (i = 1; i < 50; ip++){
		sm = 0.0;
		for (ip = 0; ip < M_Size - 1; ip++){
			for (iq = ip + 1; iq < M_Size; iq++)
				sm += abs(a[ip][ip]);
		}
		if (sm == 0.0){
			eigsrt(d, v);
			return 0;
		}
		else
			tresh = 0.0;
		for (ip = 0; ip < M_Size - 1; ip++){
			for (iq = ip + 1; iq < M_Size; iq++){
				g = 100.0*abs(a[ip][iq]);
				if (i>4 && g <= EPS*abs(d[ip]) && g <= EPS*abs(d[iq]))
					a[ip][iq] = 0.0;
				else if (abs(a[ip][iq])>tresh){
					h = d[iq] - d[ip];
					if (g <= EPS*abs(h))
						t = (a[ip][iq]) / h;
					else {
						theta = 0.5*h / (a[ip][iq]);
						t = 1.0 / (abs(theta) + sqrt(1.0 + sqr(theta)));
						if (theta < 0.0)t = -t;
					}
					c = 1.0 / sqrt(1 + sqr(t));
					s = t*c;
					tau = s / (1.0 + c);
					h = t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
					for (j = 0; j < ip; j++)
						rot(a, s, tau, j, ip, j, iq);
					for (j = ip + 1; j < iq; j++)
						rot(a, s, tau, ip, j, j, iq);
					for (j = iq + i; j < M_Size; j++)
						rot(a, s, tau, ip, j, iq, j);
					for (j = 0; j < M_Size; j++)
						rot(v, s, tau, j, ip, j, iq);
					++nrot;
				}
			}
		}
		for (ip = 0; ip < M_Size; ip++){
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	return nrot;
}


void rot(double(**a), const double s, const double tau, const int i, const int j, const int k, const int l){
	double g = a[i][j];
	double h = a[k][l];
	a[i][j] = g - s*(h + g*tau);
	a[k][l] = h + s*(g - h*tau);
}
void eigsrt(double(*d), double(**v)){
	int k;
	int n = M_Size;
	for (int i = 0; i < n - 1; i++){
		double p = d[k = i];
		for (int j = i; j < n; j++)
		if (d[j] >= p)p = d[k = j];
		if (k != i){
			d[k] = d[i];
			d[i] = p;
			if (v != NULL)
			for (int j = 0; j < n; j++){
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}

