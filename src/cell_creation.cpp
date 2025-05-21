#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

bool triangulation(vector<unsigned int>& input_vertex,vector<double>& verteces)
{
	
}

vector<double> find_vertex(vector<unsigned int>& v)
{
	
	unsigned int p=v[0];
	unsigned int q=v[1];
	unsigned int b=v[2];
	unsigned int c=v[3];
	
	if (p!=3)
		return false;
	
	if (q<3 or q>5)
		return false;
	
	if( q==3)
	{
		std::vector<std::vector<double>> tetra_vert={{1.0,1.0,1.0},{-1.0,-1.0,1.0},{-1.0,1.0,-1.0},{1.0,-1.0,-1.0}};
		return tetra_vert;
	};
	
	if( q==4)
	{
		std::vector<std::vector<double>> octa_vert={{1.0,0.0,0.0},{-1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,-1.0,0.0},{0.0,0.0,1.0},{0.0,0.0,-1.0}};
		return octa_vert;
	};
	
	if( q==5)
	{
		double golden_ratio = (1+sqrt(5))/2;
		std::vector<std::vector<double>> ico_vert={{1.0,golden_ratio,0.0},{-1.0,golden_ratio,0.0},{1.0,- golden_ratio,0.0},{-1.0,- golden_ratio,0.0},{golden_ratio,0.0,1.0},{-golden_ratio,0.0,1.0},{golden_ratio,0.0,-1.0},{-golden_ratio,0.0,-1.0},{0.0,1.0,golden_ratio},{0.0,-1.0,golden_ratio},{0.0,1.0,-golden_ratio},{0.0,-1.0,-golden_ratio}};
		return ico_vert;
	};
	
}