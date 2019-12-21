#include "pch.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include "Gaus.h"

using namespace std;

void Nev(vector<double>&, vector<double>, vector<double>, double a, double kk, double);
void Jac(vector< vector<double>>&, vector<double>, double a, double kk, double);
void Jac2(vector< vector<double>>&, vector<double>, double);
double fd1(vector<double> x, vector<double> y, double a, double kk, double h);
double fd2(vector<double>, vector<double>);

bool Hiuton(double a, double kk, vector<double>& x, vector<double> y, double h)
{
	int n = 3;
	vector<double> x1;
	vector<double> deltx;
	vector<double> F;
	vector<int> p;
	vector<vector<double>> J;
	J.resize(n);
	for (int i = 0; i < n; i++)
	{
		F.push_back(0);
		x1.push_back(0);
		x1.push_back(0);
		deltx.push_back(0);
		p.push_back(i);
		for (int j = 0; j < n; j++)
			J[i].push_back(0);
	}
	
	double e1 = 0.000000001, e2 = 0.000000001, d1 = 1, d2 = 1;
	int k = 1, NIT = 100;
	//cout << "k				d1				d2" << endl;


	while ((d1 > e1 || d2 > e2) && k < NIT)
	{
		Nev(F, x, y, a, kk, h);
		Jac(J, x, a, kk, h);
		//Jac2(J, x, 0.01);
		if (!Gaus(p, deltx, J, F))
			return 0;
		for (int i = 0; i < n; i++)
			x1[i] = x[i] + deltx[i];
		d1 = fd1(x,y,a,kk,h);
		d2 = fd2(x, x1);

		//cout << setprecision(20) << fixed << k << "				" << d1 << "		" << d2 << endl;
		k++;
		x = x1;
		if (k >= NIT)
			cout << "ERROR: IER = 2" << endl;
	}

	//Out(x);

}

void Nev(vector<double>& F, vector<double> x, vector<double> y, double a, double kk, double h)
{
	F[0] = x[0] - y[0] - h * ((kk - a)*x[1] * x[2]) / a;
	F[1] = x[1] - y[1] - h * ((kk + a)*x[0] * x[2]) / kk;
	F[2] = x[2] - y[2] - h * ((a - kk)*x[0] * x[1]) / a;
}

void Jac(vector< vector<double>>& J, vector<double> x, double a, double kk, double h)
{
	J[0][0] = 1;
	J[0][1] = (h*(kk - a)*x[2]) / a;
	J[0][2] = (h*(kk - a)*x[1]) / a;
	J[1][0] = (h*(kk + a)*x[2]) / kk;
	J[1][1] = 1;
	J[1][2] = (h*(kk + a)*x[0]) / kk;
	J[2][0] = (h*(kk + a)*x[1]) / kk;
	J[2][1] = (h*(kk + a)*x[0]) / kk;
	J[2][2] = 1;
}

void Jac2(vector< vector<double>>& J, vector<double> x, double M)
{
	J[0][0] = ((1.5 * pow(x[0] + M * x[0], 3) - pow(x[1], 2) - 1) - (1.5 * pow(x[0], 3) - pow(x[1], 2) - 1)) / (M * x[0]);
	J[0][1] = ((1.5 * pow(x[0], 3) - pow(x[1] + M * x[1], 2) - 1) - (1.5 * pow(x[0], 3) - pow(x[1], 2) - 1)) / (M * x[1]);
	J[1][0] = (((x[0] + M * x[0]) * pow(x[1], 3) - x[1] - 4) - (x[0] * pow(x[1], 3) - x[1] - 4)) / (M * x[0]);
	J[1][1] = ((x[0] * pow(x[1] + M * x[1], 3) - (x[1] + M * x[1]) - 4) - (x[0] * pow(x[1], 3) - x[1] - 4)) / (M * x[1]);
}

double fd1(vector<double> x, vector<double> y, double a, double kk, double h)
{
	double x1 = x[0] - y[0] - h * ((kk - a)*x[1] * x[2]) / a;
	double x2 = x[1] - y[1] - h * ((kk + a)*x[0] * x[2]) / kk;
	double x3 = x[2] - y[2] - h * ((a - kk)*x[0] * x[1]) / a;

	double max = 0;
	if (abs(x1) > abs(x2))
		max = abs(x1);
	else
		max = abs(x2);
	if (abs(x3) > max)
		max = abs(x3);

	return max;
}

double fd2(vector<double> x, vector<double> x1)
{
	int n = x.size();
	double max = 0;
	for (int i = 0; i < n; i++)
	{
		if (abs(x1[i]) < 1)
		{
			if (abs(x1[i] - x[i]) > max)
				max = abs(x1[i] - x[i]);
		}
		else
		{
			if (abs((x1[i] - x[i]) / x1[i]) > max)
				max = abs((x1[i] - x[i]) / x1[i]);
		}
	}

	return max;
}


