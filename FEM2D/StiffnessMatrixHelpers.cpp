#include <vector>
#include "HelperFunctions.h"
#include "iostream"

using namespace std;

void LocalStiffnessMatrixHelper(const double& s, const double& n, const int& i, const vector<int>& elem, const vector<double>& node, double* Qe, double& detJac)
{
	double *N, *dNds, *dNdn;
	N = new double[4];
	dNds = new double[4];
	dNdn = new double[4];
	ShapeFunctionGenerator(s, n, N, dNds, dNdn);
	double x[4], y[4];
	int temp;
	for (int j = 0; j < 4; j++)
	{
		temp = 2 * elem[4 * i + j];
		x[j] = node[2 * elem[4 * i + j]];
		y[j] = node[2 * elem[4 * i + j] + 1];
	}

	double Jac[4] = { 0 };
	for (int j = 0; j < 4; j++)
	{
		Jac[0] += dNds[j] * x[j];
		Jac[1] += dNds[j] * y[j];
		Jac[2] += dNdn[j] * x[j];
		Jac[3] += dNdn[j] * y[j];
	}

	detJac = Jac[0] * Jac[3] - Jac[1] * Jac[2];
	double invJac[4];
	invJac[0] = 1.0 / detJac * Jac[3];
	invJac[1] = -1.0 / detJac * Jac[1];
	invJac[2] = -1.0 / detJac * Jac[2];
	invJac[3] = 1.0 / detJac * Jac[0];

	for (int j = 0; j < 24; j++)
		Qe[j] = 0;

	double dNdx[4], dNdy[4];
	for (int j = 0; j < 4; j++)
	{
		dNdx[j] = invJac[0] * dNds[j] + invJac[1] * dNdn[j];
		dNdy[j] = invJac[2] * dNds[j] + invJac[3] * dNdn[j];
		Qe[2 * j] = dNdx[j];
		Qe[16 + 2 * j + 1] = dNdx[j];
		Qe[8 + 2 * j + 1] = dNdy[j];
		Qe[16 + 2 * j] = dNdy[j];
	}

	delete[] N;
	delete[] dNds;
	delete[] dNdn;
}

void ShapeFunctionGenerator(const double& s, const double& n, double* N, double* dNds, double* dNdn)
{
	N[0] = 1.0 / 4 * (1 - s)*(1 - n);
	N[1] = 1.0 / 4 * (1 + s)*(1 - n);
	N[2] = 1.0 / 4 * (1 + s)*(1 + n);
	N[3] = 1.0 / 4 * (1 - s)*(1 + n);

	dNds[0] = -1.0 / 4 * (1 - n);
	dNds[1] = 1.0 / 4 * (1 - n);
	dNds[2] = 1.0 / 4 * (1 + n);
	dNds[3] = -1.0 / 4 * (1 + n);

	dNdn[0] = -1.0 / 4 * (1 - s);
	dNdn[1] = -1.0 / 4 * (1 + s);
	dNdn[2] = 1.0 / 4 * (1 + s);
	dNdn[3] = 1.0 / 4 * (1 - s);
}

void ComputeLocalStiffnesMatrix(double* K_Local, const double* Qe, const double* Weight, const int& j1, const int& j2, const double* cons_matrix, const double& detJac, double& thickness)
{
	double intermediate1[24];
	for (int i = 0; i < 8; i++)
	{
		intermediate1[i] = *(cons_matrix + 0) * Qe[i] + *(cons_matrix + 1) * Qe[8 + i] + *(cons_matrix + 2) * Qe[16 + i];
		intermediate1[8 + i] = *(cons_matrix + 3) * Qe[i] + *(cons_matrix + 4) * Qe[8 + i] + *(cons_matrix + 5) * Qe[16 + i];
		intermediate1[16 + i] = *(cons_matrix + 6) * Qe[i] + *(cons_matrix + 7) * Qe[8 + i] + *(cons_matrix + 8) * Qe[16 + i];
	}

	double intermediate2[64];
	int i1, i2;
	for (int i = 0; i < 64; i++)
	{
		i1 = i / 8;
		i2 = i % 8;
		intermediate2[i] = Qe[i1] * intermediate1[i2] + Qe[i1 + 8] * intermediate1[i2 + 8] + Qe[i1 + 16] * intermediate1[i2 + 16];
		K_Local[i] += Weight[j1] * Weight[j2] * intermediate2[i] * thickness * detJac;
	}
	
}

void GlobalMatrixAllocator(double *K_Local, vector<double> &K, vector<int>& elem, const int& num_nodes, const int& i)
{
	int jL, jG, kL, nodej, nodek;
	for (int j = 0; j < 64; j++)
	{
		jL = j / 8;
		kL = j % 8;
		nodej = elem[4 * i + (jL / 2)];
		nodek = elem[4 * i + (kL / 2)];
		jG = (2 * nodej + (jL % 2)) * 2 * num_nodes + 2 * nodek + (kL % 2);
		K[jG] += K_Local[j];
	}
}

void ShapeFunctionSurface(double& s, double *Nt, double *dNtds)
{
	Nt[0] = (1 - s) / 2;
	Nt[1] = (1 + s) / 2;

	dNtds[0] = -1.0 / 2;
	dNtds[1] = 1.0 / 2;
}