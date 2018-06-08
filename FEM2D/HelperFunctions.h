#pragma once
#ifndef HELPERFUNCTION_H
#define HELPERFUNCTION_H

#include <vector>
using namespace std;

void LocalStiffnessMatrixHelper(const double& s, const double& n, const int& i, const vector<int>& elem, const vector<double>& node, double* Qe, double& detJac);
void ShapeFunctionGenerator(const double& s, const double& n, double* N, double* dNds, double* dNdn);
void ComputeLocalStiffnesMatrix(double* K_Local, const double* Qe, const double* Weight, const int& j1, const int& j2, const double* cons_matrix, const double& detJac, double& thickness);
void GlobalMatrixAllocator(double *K_Local, vector<double> &K, vector<int>& elem, const int& num_nodes, const int& i);
void ShapeFunctionSurface(double& s, double *Nt, double *dNtds);

#endif // !HELPERFUNCTION_H