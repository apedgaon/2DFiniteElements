#pragma once
#ifndef PRINTHEADER_H
#define	PRINTHEADER_H

#include <fstream>
#include <vector>
#include <cstring>
using namespace std;

template <typename T>
void printv(vector<T> const &v, string& wsFileName)
{
	ofstream write(wsFileName);
	for (int i = 0; i < v.size(); i++)
		write << v[i] << ",";
	write.close();
}

template <class T>
void printv(const T &v, const int &size)
{
	ofstream write("temp.dat");
	for (int i = 0; i < size; i++)
		write << v[i] << ",";
	write.close();
}

template <typename T>
void printm(vector<T> const &v, const int& col)
{
	ofstream write("temp.dat");
	for (int i = 0; i < v.size(); i++)
	{
		if ((i + 1) % col != 0)
			write << v[i] << ",";
		else
			write << v[i] << "\n";
	}
	write.close();
}

template <typename T>
void printm(const T* v, const int row, const int col)
{
	ofstream write("temp.dat");
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			write << v[col*i + j] << ",";
		}
		write << "\n";
	}
}

#endif