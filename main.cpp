#include <algorithm>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <mutex>
#include <netcdf>
#include <numeric>
#include <vector>
#include <filesystem>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnUserParameters.h"

#include "fitting/projectfit.hpp"

#include "grapher.hpp"
#include "oQTM/oQTM.hpp"
#include "speed_of_sound.hpp"

namespace fs = std::filesystem;

typedef std::vector<float> argodata_t;
typedef std::vector<argodata_t> argotable_t;

struct argodata_struct {
  std::size_t N_PROF;
  argotable_t depths;
  argotable_t speeds_of_sound;
  argodata_t lats;
  argodata_t longs;
  argodata_t dates;
  std::vector<std::string> platforms;
};

const argodata_struct read_file_data(fs::path filepath) {
  using namespace netCDF;
  NcFile datafile(filepath.string(), NcFile::read);

  static const std::size_t N_PROF = datafile.getDim("N_PROF").getSize();
  static const std::size_t N_LEVELS = datafile.getDim("N_LEVELS").getSize();

  argodata_struct data{0,
                       argotable_t(),
                       argotable_t(),
                       argodata_t(N_PROF),
                       argodata_t(N_PROF),
                       argodata_t(N_PROF),
                       std::vector<std::string>()};

  data.depths.reserve(N_PROF);
  data.speeds_of_sound.reserve(N_PROF);
  data.platforms.reserve(N_PROF);

  argodata_t lats(N_PROF);
  argodata_t longs(N_PROF);
  argodata_t dates(N_PROF);
  std::vector<std::string> platforms;

  NcVar latsVar, longsVar, dateVar, platformVar;
  latsVar = datafile.getVar("LATITUDE");
  longsVar = datafile.getVar("LONGITUDE");
  dateVar = datafile.getVar("JULD");
  platformVar = datafile.getVar("PLATFORM_NUMBER");

  std::vector<std::array<char, 8>> names(N_PROF);
  platformVar.getVar(names.data());

  for (auto name : names) {
    auto strnm = std::string(name.begin(), name.end());
    std::string::iterator end_pos =
        std::remove(strnm.begin(), strnm.end(), ' ');
    platforms.push_back(std::string(strnm.begin(), end_pos));
  }

  NcVar tempVar, depthVar, salinityVar;
  tempVar = datafile.getVar("TEMP");
  depthVar = datafile.getVar("DEPH_CORRECTED");
  salinityVar = datafile.getVar("PSAL_CORRECTED");

  if (tempVar.isNull() || depthVar.isNull() || salinityVar.isNull())
    exit(1);

  latsVar.getVar(lats.data());
  longsVar.getVar(longs.data());
  dateVar.getVar(dates.data());

  std::mutex data_mutex;

  for (std::size_t j = 0; j < N_PROF; j++) {
    auto tempIn = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
    auto salIn = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
    auto depthIn = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
    auto speed_of_sound_vec =
        std::vector<float_t>(N_LEVELS, std::nan("nodata"));

    std::vector<std::size_t> start{j, 0};
    std::vector<std::size_t> count{1, N_LEVELS};

    salinityVar.getVar(start, count, salIn.data());
    tempVar.getVar(start, count, tempIn.data());
    depthVar.getVar(start, count, depthIn.data());

#pragma omp parallel for
    for (std::size_t i = 0; i < N_LEVELS; i++) {
      if (tempIn[i] != 99999 && salIn[i] != 99999 && depthIn[i] != 99999) {
        speed_of_sound_vec[i] = speed_of_sound::speed_of_sound(
            speed_of_sound::pressure_at_depth(depthIn[i], data.lats[j]),
            tempIn[i], salIn[i]);
      } else {
        tempIn[i] = std::nan("nodata");
        salIn[i] = std::nan("nodata");
        depthIn[i] = std::nan("nodata");
      }
    }

    {
      auto it = std::begin(speed_of_sound_vec);
      for (; it != std::end(speed_of_sound_vec); ++it) {
        if (std::isnan(*it))
          break;
      }
      auto size =
          (long unsigned int)std::distance(std::begin(speed_of_sound_vec), it);

      if (size == 0) {
        continue;
      }
      speed_of_sound_vec.resize(size);
      depthIn.resize(size);
      {
        std::lock_guard<std::mutex> guard(data_mutex);
        data.lats.push_back(lats[j]);
        data.longs.push_back(longs[j]);
        data.dates.push_back(dates[j]);
        data.platforms.push_back(platforms[j]);
        data.N_PROF++;

        data.depths.push_back(depthIn);
        data.speeds_of_sound.push_back(speed_of_sound_vec);
      }
    }
  }
  return std::move(data);
}

int other_main(int argc, char *argv[]) {

  std::cout << "ARGO (Met Office Hadley Centre) Data" << std::endl;

  if (argc == 1) {
    std::cout << "Please pass a filename to the programme" << std::endl;
    return 1;
  }

  fs::path argument(argv[1]);

  auto data = read_file_data(argument);

  std::cout << "[INFO]: FINISHED PROCESSING\n" << std::endl;
#pragma omp parallel for
  for (std::size_t i = 0; i < data.N_PROF; i++) {
    auto sos_vect = data.speeds_of_sound[i];
    auto depth_vect = data.depths[i];

    std::vector<double> doublevect(sos_vect.begin(), sos_vect.end());

    Chisquared fcn(std::vector<double>(depth_vect.begin(), depth_vect.end()),
                   doublevect);
    ROOT::Minuit2::MnUserParameters upar;
    // upar.Add("x", (double) 500, 1);
    upar.Add("y", (double)doublevect[0], 1);
    upar.Add("m_1", 0, 1);
    upar.Add("m_2", 0, 1);

    ROOT::Minuit2::MnMigrad migrad(fcn, upar);
    ROOT::Minuit2::FunctionMinimum min = migrad();

    std::vector<double> fitted_to_params = fcn.fitted_to_minimisation(min);
    std::string filename = "images/" + data.platforms[i] + ".eps";
    Gnuplot gp;
    gp << "set terminal postscript enhanced eps color size 8,3\n"
          "set output '"
       << filename << "'\n";
    gp << "set multiplot layout 2,1\n"
       << "set key off\n";
    grapher::plot_lines(gp, depth_vect, doublevect, sos_vect, fitted_to_params);
    grapher::plot_lines(gp, moving_average(depth_vect),
                        moving_average(sos_vect));
  }
  return 0;
}

template <typename MESH>
void read_to_tree(fs::path filepath, MESH &mesh){
  using namespace netCDF;

  NcFile datafile(filepath.string(), NcFile::read);

  static const std::size_t N_PROF = datafile.getDim("N_PROF").getSize();
  static const std::size_t N_LEVELS = datafile.getDim("N_LEVELS").getSize();



}

int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "Please pass a filename to the programme" << std::endl;
    return 1;
  }

  fs::path argument(argv[1]);

  auto data = read_file_data(argument);

  oQTM_Mesh<float_t, 3> globemesh;
  oQTM_Mesh<float_t, 3>::location_t l{0, 0, 0};
  oQTM_Mesh<float_t, 3>::location_t l_1{0, 0, 3};

  globemesh.insert(l, 0.15, 2.0);
  globemesh.insert(l_1, 45, 3.0);

  auto a = globemesh.get_points(std::vector<size_t>{0, 0, 0});
  for (auto i : a) {
    std::cout << i.first << " " << i.second << std::endl;
  }
}