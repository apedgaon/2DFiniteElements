#include <fstream>
#include <cstring>
#include <vector>
#include <cassert>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])
{
	string path = "D:/VSProj/2DFiniteElements/FEM2D/";	
	ifstream read_file;
	ofstream write;

	// Nodal Coordinates Reading
	string pathCoord = path;
	pathCoord.append("Coord.dat");	
	vector<double> node;
	double coordinates;
	read_file.open(pathCoord);
	assert(read_file.is_open());
	while (!read_file.eof())
	{
		read_file >> coordinates;
		node.push_back(coordinates);
	}

	read_file.close();
	read_file.clear();
	int num_nodes = node.size() / 2;

	// Connectivity Reading
	string  pathConnectivity;
	pathConnectivity = path;
	pathConnectivity.append("Connectivity.dat");
	vector<int> elem;
	int elemID;	
	read_file.open(pathConnectivity);
	assert(read_file.is_open());
	while (!read_file.eof())
	{
		read_file >> elemID;
		elem.push_back(elemID - 1);
	}

	read_file.close();
	read_file.clear();
	int num_elem = elem.size() / 4;

	// Writing "DBC.dat"
	string pathDir = path;
	pathDir.append("DBC.dat");
	write.open(pathDir);
	assert(write.is_open());
	double ux = 0;
	double uy = 0;
	for (int i = 0; i < num_nodes; i++)
	{
		if (fabs(node[2 * i] - 0) < 1e-5)
		{
			write << i << "\t" << 0 << "\t" << ux << "\n";			
		}

		if (fabs(node[2 * i + 1] - 0) < 1e-5)
		{
			write << i << "\t" << 1 << "\t" << uy << "\n";
		}
	}
	write.close();
	write.clear();

	// Wrinting "NBC.dat"
	string pathNeum = path;
	pathNeum.append("NBC.dat");
	write.open(pathNeum);
	assert(write.is_open());
	int quad_edge[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };
	double tx = 0.5;
	for (int i = 0; i < num_elem; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (fabs(node[2 * elem[4 * i + quad_edge[2 * j]]] - 2) < 1e-5 && fabs(node[2 * elem[4 * i + quad_edge[2 * j + 1]]] - 2) < 1e-5)
			{
				write << i << "\t" << quad_edge[2 * j] << "\t" << quad_edge[2 * j + 1] << "\t" << 0 << "\t" << tx << "\n";
			}
		}
	}

	write.close();
	write.clear();

	return 0;
}