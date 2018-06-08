#include <vector>
#include <cmath>
#include <iostream>
#include "Solver.h"

using namespace std;

void CompressedRowStorage(vector<double>& A, vector<double>& vdCompressedVector, vector<int>& viRowPointer, vector<int>& viColumnIndex, const int& num_nodes)
{
	int counter = 0;
	int column;
	for (int i = 0; i < num_nodes*num_nodes; i++)
	{
		column = i % num_nodes;
		if (column == 0)
			viRowPointer.push_back(counter);

		if (fabs(A[i] - 0.0) > 1e-5)
		{
			vdCompressedVector.push_back(A[i]);
			viColumnIndex.push_back(column);
			counter++;
		}
	}
	viRowPointer.push_back(counter);
}

void CGSolver(vector<double>& A, vector<double>& u, vector<double>& b, const int& n)
{
	vector<double> A_Compressed;
	vector<int> row_pointer;
	vector<int> col_index;
	CompressedRowStorage(A, A_Compressed, row_pointer, col_index, n);

	vector<double> dResidual(n), dPreconditioner(n), dPreCondResidual(n);
	dResidual = b;
	double dNormResidual = 0;
	double dNormb = 0;
	for (int i = 0; i < n; i++)
	{
		dNormResidual += pow(dResidual[i], 2);
		dNormb += pow(b[i], 2);
		dPreconditioner[i] = A[i*(n + 1)];
		dPreCondResidual[i] = 1.0 / dPreconditioner[i] * dResidual[i];
	}

	vector<double> dConjGrad(n);
	dConjGrad = dPreCondResidual;
	int k = 0;
	double alpha, alpha_numerator, alpha_denominator, beta, beta_numerator, beta_denominator;
	vector<double> alpha_denominator_sub(n);
	while (k < 1e5)
	{
		alpha_numerator = 0;
		alpha_denominator = 0;
		for (int i = 0; i < n; i++)
		{
			alpha_numerator += dResidual[i] * dPreCondResidual[i];
			alpha_denominator_sub[i] = 0;
			for (int j = row_pointer[i]; j < row_pointer[i + 1]; j++)
				alpha_denominator_sub[i] += A_Compressed[j] * dConjGrad[col_index[j]];

			alpha_denominator += dConjGrad[i] * alpha_denominator_sub[i];
		}

		alpha = alpha_numerator / alpha_denominator;
		beta_denominator = alpha_numerator;
		for (int i = 0; i < n; i++)
		{
			u[i] += alpha*dConjGrad[i];
			dResidual[i] -= alpha*alpha_denominator_sub[i];
		}

		dNormResidual = 0;
		for (int i = 0; i < n; i++)
			dNormResidual += pow(dResidual[i], 2);

		if (dNormResidual < 1e-5 * dNormb)
			break;

		beta_numerator = 0;
		for (int i = 0; i < n; i++)
		{
			dPreCondResidual[i] = 1.0 / dPreconditioner[i] * dResidual[i];
			beta_numerator += dPreCondResidual[i] * dResidual[i];
		}

		beta = beta_numerator / beta_denominator;

		for (int i = 0; i < n; i++)
			dConjGrad[i] = dPreCondResidual[i] + beta*dConjGrad[i];

		k++;
	}

	if (k == 1e5)
		cerr << "******* Convergence did not occur !" << "\n";

	cout << k << "\n";
	getchar();
}