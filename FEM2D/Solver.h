#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
using namespace std;

void CompressedRowStorage(vector<double>&, vector<double>&, vector<int>&, vector<int>&, const int&);
void CGSolver(vector<double>&, vector<double>&, vector<double>&, const int&);

#endif // !SOLVER_H
