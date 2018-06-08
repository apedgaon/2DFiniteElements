#include <iostream>
#include <vector>
#include <cmath>
#include "BoundaryConditions.h"
#include "HelperFunctions.h"

using namespace std;

Dirichlet::Dirichlet()
{
	ndID = 0;
	dof = 0;
	val = 0;
	size = 0;
}

Dirichlet::Dirichlet(const int& n)
{
	size = n;
	ndID = new int[size];
	dof = new int[size];
	val = new double[size];
}

Dirichlet::~Dirichlet()
{
	if (ndID != 0)
	{
		delete[] ndID;
		delete[] dof;
		delete[] val;
	}
}

void Dirichlet::inputDir(const int& a, const int& b, const double& c)
{
	if (size != 0)
	{
		Dirichlet temp(size);
		for (int i = 0; i < size; i++)
		{
			temp.ndID[i] = ndID[i];
			temp.dof[i] = dof[i];
			temp.val[i] = val[i];
		}
		delete[] ndID;
		delete[] dof;
		delete[] val;

		size++;
		ndID = new int[size];
		dof = new int[size];
		val = new double[size];
		for (int i = 0; i < size - 1; i++)
		{
			ndID[i] = temp.ndID[i];
			dof[i] = temp.dof[i];
			val[i] = temp.val[i];
		}
		
		ndID[size - 1] = a;
		dof[size - 1] = b;
		val[size - 1] = c;
	}
	else
	{
		size++;
		ndID = new int[size];
		dof = new int[size];
		val = new double[size];
		ndID[size - 1] = a;
		dof[size - 1] = b;
		val[size - 1] = c;
	}
}

int Dirichlet::getsize()
{
	return size;
}

void Dirichlet::print()
{
	for (int i = 0; i < size; i++)
		std::cout << ndID[i] << "," << dof[i] << "," << val[i] << "\n";

	getchar();
}

void Dirichlet::ApplyDirichlet(vector<double>& K, vector<double>& F, const int& num_nodes)
{
	int iG, jG, kG;
	for (int i = 0; i < size; i++)
	{
		jG = 2 * ndID[i] + dof[i];
		iG = jG * 2 * num_nodes;
		for (int j = 0; j < 2 * num_nodes; j++)
			K[iG + j] = 0;

		K[iG + jG] = 1;
		F[jG] = val[i];
		for (int k = 0; k < 2 * num_nodes; k++)
		{
			kG = k * 2 * num_nodes + jG;
			if (kG != iG + jG)
			{
				F[k] -= K[kG] * F[jG];
				K[kG] = 0;
			}
		}
	}
}


Neumann::Neumann()
{
	elemID = 0;
	node1 = 0;
	node2 = 0;
	dof = 0;
	val = 0;
	size = 0;
}

Neumann::Neumann(const int& n)
{
	size = n;
	elemID = new int[size];
	node1 = new int[size];
	node2 = new int[size];
	dof = new int[size];
	val = new double[size];
}

Neumann::~Neumann()
{
	delete[] elemID;
	delete[] node1;
	delete[] node2;
	delete[] dof;
	delete[] val;
	size = 0;
}

void Neumann::inputNeumann(const int& elID, const int& nd1, const int& nd2, const int& df, double& value)
{
	if (size != 0)
	{
		Neumann temp(size);
		for (int i = 0; i < size; i++)
		{
			temp.elemID[i] = elemID[i];
			temp.node1[i] = node1[i];
			temp.node2[i] = node2[i];
			temp.dof[i] = dof[i];
			temp.val[i] = val[i];
		}

		
		delete[] elemID;
		delete[] node1;
		delete[] node2;
		delete[] dof;
		delete[] val;
		//*this = temp;/*
		size++;
		elemID = new int[size];
		node1 = new int[size];
		node2 = new int[size];
		dof = new int[size];
		val = new double[size];
		
		
		for (int i = 0; i < size - 1; i++)
		{
			elemID[i] = temp.elemID[i];
			node1[i] = temp.node1[i];
			node2[i] = temp.node2[i];
			dof[i] = temp.dof[i];
			val[i] = temp.val[i];
		}

		elemID[size - 1] = elID;
		node1[size - 1] = nd1;
		node2[size - 1] = nd2;
		dof[size - 1] = df;
		val[size - 1] = value;
	}
	else
	{
		size++;
		elemID = new int[size];
		node1 = new int[size];
		node2 = new int[size];
		dof = new int[size];
		val = new double[size];
		elemID[size - 1] = elID;
		node1[size - 1] = nd1;
		node2[size - 1] = nd2;
		dof[size - 1] = df;
		val[size - 1] = value;
	}
}

int Neumann::getsize()
{
	return size;
}

void Neumann::print()
{
	for (int i = 0; i < size; i++)
		std::cout << elemID[i] << "," << node1[i] << "," << node2[i] << "," << dof[i] << "," << val[i] << "\n";

	getchar();
}

void Neumann::ApplyNeum(vector<int>& elem, vector<double>& node, vector<double>& F, double& thickness)
{
	int jG, kG;
	double JacS;
	double s[2] = { -0.577350269189626, 0.577350269189626 };
	double Nt[2], dNtds[2];
	double Qt[8], temp[16], F_Local[4], traction[4];
	int k1, k2;
	for (int i = 0; i < size; i++)
	{
		jG = elem[4 * elemID[i] + node1[i]];
		kG = elem[4 * elemID[i] + node2[i]];
		JacS = sqrt((pow((node[2 * jG] - node[2 * kG]), 2) + pow((node[2 * jG + 1] - node[2 * kG + 1]), 2))) / 2;
		
		for (int j = 0; j < 16; j++)
			temp[j] = 0;

		for (int j = 0; j < 2; j++)
		{
			ShapeFunctionSurface(s[j], Nt, dNtds);
			for (int k = 0; k < 8; k++)
				Qt[k] = 0;

			Qt[0] = Nt[0];
			Qt[2] = Nt[1];
			Qt[5] = Nt[0];
			Qt[7] = Nt[1];
			for (int k = 0; k < 16; k++)
			{
				k1 = k / 4;
				k2 = k % 4;
				temp[k] += Qt[k1] * Qt[k2] + Qt[k1 + 4] * Qt[4 + k2];
			}
		}
		
		for (int j = 0; j < 4; j++)
			traction[j] = 0;

		traction[dof[i]] = val[i];
		traction[dof[i] + 2] = val[i];

		for (int j = 0; j < 4; j++)
		{
			F_Local[j] = thickness * JacS * (temp[4 * j] * traction[0] + temp[4 * j + 1] * traction[1] + temp[4 * j + 2] * traction[2] + temp[4 * j + 3] * traction[3]);
		}

		F[2 * jG] += F_Local[0];
		F[2 * jG + 1] += F_Local[1];
		F[2 * kG] += F_Local[2];
		F[2 * kG + 1] += F_Local[3];
	}
}