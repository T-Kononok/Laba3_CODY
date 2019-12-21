#include "pch.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include "Hiuton.h"

using namespace std;

int a = 1, k = 2;

void uf(vector<double>& uu, vector<double>& u);
void uh(vector<double>& uu, vector<double>& u, double h);
double shag(double e, double hmax, vector<double>&);
void ekf(vector<double>& ek, vector<double> uk_1, vector<double> u, vector<double> uk1, double h, double hk_1);
void u_uu(vector<double>& u, vector<double>& uu);
double shag1(double h, double e, vector<double> ek, bool);

int main()
{
	setlocale(LC_ALL, "Russian");
	bool v;

	cout << "Явный - 1, неявный - 0 : ";
	cin >> v;

	double t = 0, T = 1, e = pow(10, -3);
	vector<double> u;
	u.assign(3, 1);

	if (v)
	{
		double hmax = 1, h;
		vector<double> uu;
		uu.assign(3, 0);
		int i = 0;
		do
		{
			i++;
			uf(uu, u);
			h = shag(e, hmax, uu);
			uh(uu, u, h);
			t += h;
			cout << setprecision(6) << fixed << i << " " << u[0] << " " << u[1] << " " << u[2] << " " << t << endl;
		} while (t < T);
	}
	else
	{
		double hmax = 1, hmin = 0.01, h = hmin, hk1 = hmin, hk_1 = hmin, tk1;
		vector<double> uu, uk1, uk_1, ek;
		uk1.assign(3, 0);
		uk_1.assign(3, 0);
		ek.assign(3, 0);
		uu.assign(3, 0);

		int i = 0, z = 0;
		do
		{
			do
			{
				z = 0;
				i++;
				tk1 = t + h;

				Hiuton(a, k, uk1, u, h);

				ekf(ek, uk_1, u, uk1, h, hk_1);

				if (abs(ek[0]) > e || abs(ek[1]) > e || abs(ek[2]) > e)
				{
					h /= 2;
					hk1 = h;
					u_uu(uk1,u);
					z = 1;
				}
			} while (z == 1);

			hk1 = shag1(h, e, ek, 1);/////////////////////////////////////////////////////////////

			if (hk1 > hmax)
				hk1 = hmax;

			u_uu(uk_1, u);
			u_uu(u, uk1);
			hk_1 = h;
			h = hk1;
			t = tk1;

			cout << setprecision(6) << fixed << i << " " << u[0] << " " << u[1] << " " << u[2] << " " << t << endl;
			//cout << setprecision(6) << fixed << u[2] << " ";

		} while (t < T);
	}
}

void uf(vector<double>& uu, vector<double>& u)
{
	uu[0] = ((k - a)*u[1] * u[2]) / a;
	uu[1] = ((k + a)*u[0] * u[2]) / k;
	uu[2] = ((a - k)*u[0] * u[1]) / a;
}

void uh(vector<double>& uu, vector<double>& u, double h)
{
	u[0] += h * uu[0];
	u[1] += h * uu[1];
	u[2] += h * uu[2];
}

void u_uu(vector<double>& u, vector<double>& uu)
{
	u[0] = uu[0];
	u[1] = uu[1];
	u[2] = uu[2];
}

void ekf(vector<double>& ek, vector<double> uk_1, vector<double> u, vector<double> uk1, double h, double hk_1)
{
	ek[0] = -1 * (h / (h + hk_1)) * (uk1[0] - u[0] - (h*(u[0] - uk_1[0])) / (hk_1));
	ek[1] = -1 * (h / (h + hk_1)) * (uk1[1] - u[1] - (h*(u[1] - uk_1[1])) / (hk_1));
	ek[2] = -1 * (h / (h + hk_1)) * (uk1[2] - u[2] - (h*(u[2] - uk_1[2])) / (hk_1));
}

double shag(double e, double hmax, vector<double>& uu)
{
	double h1, h2, h3, hmin;
	h1 = e / (abs(uu[0]) + (e / hmax));
	h2 = e / (abs(uu[1]) + (e / hmax));
	h3 = e / (abs(uu[2]) + (e / hmax));

	if (h1 < h2)
		hmin = h1;
	else
		hmin = h2;
	if (h3 < hmin)
		hmin = h3;

	return hmin;
}

double shag1(double h, double e, vector<double> ek, bool v)
{
	double hmin;
	vector<double> hh;
	hh.assign(3, 0);

	if (v)
	{
		for (int i = 0; i < 3; i++)
			hh[i] = pow(e / abs(ek[i]), 0.5) * h;
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			if (abs(ek[i]) > e)
				hh[i] = h / 2;
			if (e / 4 < abs(ek[i]) && abs(ek[i]) <= e)
				hh[i] = h;
			if (abs(ek[i]) < e / 4)
				hh[i] = 2 * h;
		}
	}

	if (hh[0] < hh[1])
		hmin = hh[0];
	else
		hmin = hh[1];
	if (hh[2] < hmin)
		hmin = hh[2];

	//cout << hmin << endl;
	return hmin;
}