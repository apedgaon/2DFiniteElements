#pragma once
#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <vector>
using namespace std;

class Dirichlet
{
private:
	int* ndID;
	int* dof;
	double* val;
	int size;
public:
	Dirichlet();
	Dirichlet(const int&);
	~Dirichlet();
	void inputDir(const int&, const int&, const double&);
	int getsize();
	void print();
	void ApplyDirichlet(vector<double>& K, vector<double>& F, const int& num_nodes);
};

class Neumann
{
private:
	int* elemID;
	int* node1;
	int* node2;
	int* dof;
	double* val;
	int size;
public:
	Neumann();
	Neumann(const int& n);
	~Neumann();
	void inputNeumann(const int& elID, const int& nd1, const int& nd2, const int& df, double& value);
	int getsize();
	void print();
	void ApplyNeum(vector<int>& elem, vector<double>& node, vector<double>& F, double& thickness);
};

#endif