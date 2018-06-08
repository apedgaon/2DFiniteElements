#include <fstream>
#include <vector>
#include <cassert>
#include <cstring>
#include "Print.h"
#include "BoundaryConditions.h"
#include "HelperFunctions.h"
#include "Solver.h"

using namespace std;

int main(int argc, char* argv[])
{
	// Geometric Properties
	double thickness = 0.1;
	double length = 2;
	double height = 2;

	// Material Property
	double cons_matrix[9] = { 3,1,0,1,3,0,0,0,1 };

	// Element Data Loading
	vector<int> elem;
	int elemID;
	ifstream read_file;
	read_file.open("Connectivity.dat");
	assert(read_file.is_open());
	while (!read_file.eof())
	{
		read_file >> elemID;
		elem.push_back(elemID - 1);
	}

	read_file.close();
	read_file.clear();
	int num_elem = elem.size() / 4;

	// Node Data Loading
	vector<double> node;
	double coordinates;
	read_file.open("Coord.dat");
	assert(read_file.is_open());
	while (!read_file.eof())
	{
		read_file >> coordinates;
		node.push_back(coordinates);
	}

	read_file.close();
	read_file.clear();
	int num_nodes = node.size() / 2;
	
	// Quadrature point definition
	double s[2] = { -0.577350269189626, 0.577350269189626 };
	double n[2] = { -0.577350269189626, 0.577350269189626 };
	double Weight[2] = { 1.0, 1.0 };

	// Stiffness Matrix Construction
	vector<double> K(4 * num_nodes*num_nodes,0);
	vector<double> F(2 * num_nodes,0);
	double* K_Local;
	K_Local = new double[64];
	int j1, j2;
	double *Qe, detJac;
	Qe = new double[24];
	for (int i = 1; i < num_elem; i++)
	{
		for (int j = 0; j < 64; j++)
			K_Local[j] = 0;

		for (int j = 0; j < 4; j++)
		{
			j1 = j / 2;
			j2 = j % 2;
			LocalStiffnessMatrixHelper(s[j1], n[j2], i, elem, node, Qe, detJac);
			ComputeLocalStiffnesMatrix(K_Local, Qe, Weight, j1, j2, &cons_matrix[0], detJac, thickness);
		}

		GlobalMatrixAllocator(K_Local, K, elem, num_nodes, i);
		
	}

	// Loading and Applying Neuman Boundary Conditions Loading
	Neumann NBC;
	int elID, nodeID1, nodeID2, df;
	double tvalue;
	read_file.open("NBC.dat");
	assert(read_file.is_open());
	while (!read_file.eof())
	{
		read_file >> elID >> nodeID1 >> nodeID2 >> df >> tvalue;
		NBC.inputNeumann(elID, nodeID1, nodeID2, df, tvalue);
	}
	read_file.close();
	read_file.clear();

	NBC.ApplyNeum(elem, node, F, thickness);
	
	string sFileName2("F.dat");
	printv(F, sFileName2);

	// Loading and Applying Dirichlet Boundary Conditions
	Dirichlet DBC;
	int dira, dirb;
	double dirc;
	read_file.open("DBC.dat");
	assert(read_file.is_open());
	while (!read_file.eof())
	{
		read_file >> dira >> dirb >> dirc;
		DBC.inputDir(dira, dirb, dirc);
	}
	read_file.close();
	read_file.clear();
	DBC.ApplyDirichlet(K, F, num_nodes);

	//printv(F);
	//printm(K, 2 * num_nodes);
	
	vector<double> U(2 * num_nodes, 0);
	CGSolver(K, U, F, 2 * num_nodes);
	string sFileName("U.dat");
	printv(U, sFileName);
	delete[] K_Local;
	delete[] Qe;
	return 0;
}