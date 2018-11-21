#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <netcdf>
#include <cstdio>

#include "speed_of_sound.hpp"
#include "grapher.hpp"

using std::cout;
using std::endl;



int main(int argc, char *argv[])
{
	using namespace netCDF;

	typedef std::vector<std::vector<float_t>> argodata_t;
	cout << "ARGO (Met Office Hadley Centre) Data" << endl;
	
	if (argc == 1) {
		cout << "Please pass a filename to the programme" << endl;
		return 1;
	}

	NcFile test(argv[1], NcFile::read);
	static const size_t N_PROF = test.getDim("N_PROF").getSize();
	static const size_t N_LEVELS = test.getDim("N_LEVELS").getSize();
	NcVar latsVar, longsVar, dateVar;
	latsVar = test.getVar("LATITUDE");
	longsVar = test.getVar("LONGITUDE");
	dateVar = test.getVar("JULD");

	NcVar tempVar, depthVar, salinityVar;
	tempVar = test.getVar("TEMP");
	depthVar = test.getVar("DEPH_CORRECTED");
	salinityVar = test.getVar("PSAL_CORRECTED");

	if (tempVar.isNull() || depthVar.isNull() || salinityVar.isNull())
		exit(1);

	std::vector<size_t> count{1, N_LEVELS};

	std::vector<float_t> lats(N_PROF, 0.);
	std::vector<float_t> longs(N_PROF, 0.);
	std::vector<float_t> dates(N_PROF, 0.);

	latsVar.getVar(lats.data());
	longsVar.getVar(longs.data());
	dateVar.getVar(dates.data());

	argodata_t temps(N_PROF, std::vector<float_t>(N_LEVELS, 99999));
	argodata_t salinities(N_PROF, std::vector<float_t>(N_LEVELS, 99999));
	argodata_t depths(N_PROF, std::vector<float_t>(N_LEVELS, 99999));
	argodata_t speeds_of_sound(N_PROF, std::vector<float_t>(N_LEVELS, 99999));
	for (size_t j = 0; j < N_PROF; j++)
	{
		auto tempIn = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
		auto salIn = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
		auto depthIn = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
		auto speed_of_sound_vec = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
		auto dvdx_vec = std::vector<float_t>(N_LEVELS, std::nan("nodata"));

		std::vector<size_t> start{j, 0};

		salinityVar.getVar(start, count, salIn.data());
		tempVar.getVar(start, count, tempIn.data());
		depthVar.getVar(start, count, depthIn.data());

		for (size_t i = 0; i < N_LEVELS; i++)
		{
			if (tempIn[i] != 99999 && salIn[i] != 99999 && depthIn[i] != 99999)
			{
				speed_of_sound_vec[i] = speed_of_sound::speed_of_sound(speed_of_sound::pressure_at_depth(depthIn[i], lats[j]), tempIn[i], salIn[i]);
			}
			else
			{
				tempIn[i] = std::nan("nodata");
				salIn[i] = std::nan("nodata");
				depthIn[i] = std::nan("nodata");
			}
		}
		{
			temps[j] = tempIn;
			salinities[j] = salIn;
			depths[j] = depthIn;
			speeds_of_sound[j] = speed_of_sound_vec;
		}
	}

	cout << "[INFO]: FINISHED PROCESSING\n"
			 << endl;

	for (size_t i = 0; i < N_PROF; i++)
	{
		auto sos_vect = speeds_of_sound[i];
		auto depth_vect = depths[i];

		std::vector<float_t> difference(N_LEVELS, 0);
		for (size_t j = 1; j < N_LEVELS; j++)
		{
			if (depth_vect[j] - depth_vect[j - 1] != 0 && !isnan(depth_vect[j]) && !isnan(depth_vect[j-1]))
			difference[j] = (sos_vect[j] - sos_vect[j - 1]) / (depth_vect[j] - depth_vect[j - 1]);
		}
		for ( auto diff: difference){
			cout << diff << " ";
		}
		cout << endl;
	}
	return 0;
}