#include <algorithm>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <mutex>
#include <iomanip>
#include <numeric>
#include <vector>
#include <complex>
#include <tuple>

#include "boost/filesystem.hpp"
#include <netcdf>
#include "ncException.h"
#include "ncChar.h"

#include "fitting/projectfit.hpp"

#include "grapher.hpp"
#include "oQTM/oQTM.hpp"
#include "speed_of_sound.hpp"

namespace fs = boost::filesystem;

template <typename T, typename K, typename V>
std::vector<std::tuple<T, T, K, V>> read_to_tree(const fs::path &filepath)
{
  using namespace netCDF;

  NcFile datafile(filepath.string(), NcFile::read);

  static const std::size_t N_PROF = datafile.getDim("N_PROF").getSize();
  static const std::size_t N_LEVELS = datafile.getDim("N_LEVELS").getSize();

  std::vector<std::tuple<T, T, K, V>> results;
  //results.reserve(N_PROF);

  NcVar latsVar, longsVar, dateVar, platformVar, tempVar, depthVar, salinityVar, profileQualityVar, levelsQualityVar, positionQualityVar, tempQualityVar;
  latsVar = datafile.getVar("LATITUDE");
  longsVar = datafile.getVar("LONGITUDE");
  dateVar = datafile.getVar("JULD");

  tempVar = datafile.getVar("POTM_CORRECTED");
  depthVar = datafile.getVar("DEPH_CORRECTED");
  salinityVar = datafile.getVar("PSAL_CORRECTED");

  profileQualityVar = datafile.getVar("QC_FLAGS_PROFILES");
  levelsQualityVar = datafile.getVar("QC_FLAGS_LEVELS");
  positionQualityVar = datafile.getVar("POSITION_QC");
  tempQualityVar = datafile.getVar("PROFILE_POTM_QC");

  std::vector<double_t> lats(latsVar.getDim(0).getSize(), 0), lons(longsVar.getDim(0).getSize(), 0);
  std::vector<double_t> dates(dateVar.getDim(0).getSize(), 0);

  if (tempVar.isNull() || depthVar.isNull() || salinityVar.isNull())
  {
    throw std::runtime_error("File format not as expected.");
  }
#pragma omp critical
  {
    latsVar.getVar(lats.data());
    longsVar.getVar(lons.data());
    dateVar.getVar(dates.data());
  }

  for (std::size_t i = 0; i < N_PROF; i++)
  {
    try
    {
      uint32_t profile_qc;

      char position_qc = 0;
      char potm_qc = 0;

      profileQualityVar.getVar(std::vector<size_t>{i}, &profile_qc);
      positionQualityVar.getVar(std::vector<size_t>{i}, &position_qc);
      tempQualityVar.getVar(std::vector<size_t>{i}, &potm_qc);

      if ((profile_qc & 0b11) != 0 || position_qc == 4 || potm_qc == 4)
      {
        //std::cerr << "Rejected profile " << i << " with code " << profile_qc << std::endl;
        continue;
      }
    }
    catch (const netCDF::exceptions::NcInvalidCoords &e)
    {
      std::cerr << filepath << ": " << e.what() << std::endl;
      break;
    }

    T lat = lats[i], lon = lons[i];
    K date = dates[i];

    std::vector<float_t> tempIn(tempVar.getDim(0).getSize(), std::nan("nodata"));
    std::vector<float_t> salIn(salinityVar.getDim(0).getSize(), std::nan("nodata"));
    std::vector<float_t> depthIn(depthVar.getDim(0).getSize(), std::nan("nodata"));
    std::vector<uint32_t> levelQC(levelsQualityVar.getDim(0).getSize(), std::nan("nodata"));

    std::vector<double_t> actual_depths, actual_temperatures;
    std::vector<double_t> speed_of_sound_vec;

    //actual_depths.reserve(N_LEVELS);
    //speed_of_sound_vec.reserve(N_LEVELS);
    //actual_temperatures.reserve(N_LEVELS);
    std::vector<std::size_t> start{i, 0};
    std::vector<std::size_t> count{1, tempVar.getDim(1).getSize()};

    try
    {
#pragma omp critical
      {
        levelsQualityVar.getVar(start, count, levelQC.data());
        salinityVar.getVar(start, count, salIn.data());
        tempVar.getVar(start, count, tempIn.data());
        depthVar.getVar(start, count, depthIn.data());
      }
    }
    catch (const netCDF::exceptions::NcInvalidCoords &e)
    {
      std::cerr << filepath << ": " << e.what() << std::endl;
      break;
    }

    auto no_value = [](double_t s) { return s == 99999 || std::isnan(s); };

    for (std::size_t j = 0; j < N_LEVELS; j++)
    {
      if ((levelQC[j] & 0b11) != 0)
        continue;
      auto s = speed_of_sound::speed_of_sound(
          speed_of_sound::pressure_at_depth((double_t)depthIn[j], lat),
          (double_t)tempIn[j], (double_t)salIn[j]);
      speed_of_sound_vec.push_back(s);
      actual_depths.push_back((double_t)depthIn[j]);
      actual_temperatures.push_back(tempIn[j]);
    }

    if (actual_depths.size() < 10)
    {
      continue;
    };

    speed_of_sound_vec.shrink_to_fit();
    actual_depths.shrink_to_fit();
    actual_temperatures.shrink_to_fit();

    try
    {
      auto xmin = fit::find_SOFAR_channel(speed_of_sound_vec, actual_depths);
      if (xmin > actual_depths.back()) continue;
      auto tavg = actual_temperatures[0];
#pragma omp critical
      {
        /*
        std::string filename = "images/" + std::to_string(date) + "." + std::to_string(i) + ".eps";
        Gnuplot gp;
        gp << "set terminal postscript enhanced eps color size 3,3\n"
              "set output '"
           << filename << "'\n"
           << "set key off\n";
        grapher::plot_lines(gp, actual_depths, speed_of_sound_vec); */
        const V result{xmin};

        results.push_back(std::tie(lon, lat, date, result));
      }
    }
    catch (std::runtime_error e)
    {
      continue;
    }
  }

  return results;
}

int main(int argc, char *argv[])
{
  typedef double_t value_t;
  typedef oQTM_Mesh<double_t, uint64_t, value_t, 10> mesh_t;
  /*if (argc == 1 || !fs::is_directory(argv[1]))
  {
    std::cout << "Please pass a directory or filename to the programme" << std::endl;
    return 1;
  }*/

  mesh_t globemesh;
  if (argc == 1)
    exit(1);
  auto dirs = fs::directory_iterator(argv[1]);

  std::vector<fs::directory_entry> paths(fs::begin(dirs), fs::end(dirs));
  //#pragma omp parallel
  for (std::size_t i = 0; i < paths.size(); i++)
  {
    auto path = paths[i];
    if (/*fs::is_regular_file(path) && */ path.path().extension() == ".nc")
    {
      auto points = read_to_tree<double_t, uint64_t, value_t>(path);
      //#pragma omp task
      {
        globemesh.insert(points);
      }
    }
  }

  std::cout << "[INFO]: Interpolation complete" << std::endl;

  std::vector<std::pair<double_t, double_t>> locations{
      {-40.671851, 33.792983}, // Atlantic
      {3.343785, 56.3836},     // North Sea
      {18.179810, 35.013669},  // Mediterranean
      {79.738750, -21.668352}, // Indian Ocean
      {18.179806, 35.013667},  // Bay of Bengal
      {114.008616, 15.425780}  // South China Sea
  };
    fs::create_directory("output");
    fs::create_directory("output/power_spectra");
//#pragma omp parallel for
  for (std::size_t i = 0; i < locations.size(); i++)
  {
    auto loc = globemesh.location(locations[i].first, locations[i].second);
    auto a = globemesh.get_averaged_points(loc, 4);
    std::stringstream filename;
    filename << "output/";

    for (auto i : loc)
      filename << (int)i << '.';
    filename << "results.csv";
    std::ofstream fileout(filename.str());
    fileout << "days since 1950-01-01,speed of sound minimum depth (m)" << std::endl;
    for (auto i : a)
    {
      fileout << i.first << "," << i.second << std::endl;
    }
    fileout.close();
    auto power_spectrum = fit::analyse_periodicity(a);
    filename << "output/power_spectra/";
    for (auto i : loc)
      filename << (int)i << '.';
    filename << "results.csv";
    std::ofstream fileout_ps(filename.str());
    for (auto i : power_spectrum)
    {
      fileout_ps << i.real() << "," << i.imag() << std::endl;
    }
    fileout_ps.close();
  }
}
