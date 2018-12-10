#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <netcdf>
#include <vector>
#include <fftw3.h>
#include <numeric>
#include <tuple>
#include "boost/filesystem.hpp"
#include "grapher.hpp"
#include "speed_of_sound.hpp"

namespace fs = boost::filesystem;

std::vector<float> low_pass_filter(const std::vector<float> &vector, const size_t cutoff = 15)
{
  std::vector<float> in(vector);
  float average = std::accumulate(in.begin(), in.end(), 0) / in.size();
  for (auto &i : in)
    i -= average;

  std::vector<float> out(in.size());
  fftwf_complex *fft = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * in.size() / 2 + 1);
  fftwf_plan forward = fftwf_plan_dft_r2c_1d(in.size(), in.data(), fft, FFTW_ESTIMATE);
  fftwf_execute(forward);
  fftwf_destroy_plan(forward);

  for (size_t i = 0; i < in.size() / 2 + 1; i++)
  {
    float window = 1.0 ? i < cutoff : std::sqrtf(std::expf(-(i - cutoff) * (i - cutoff) / 10.0));
    fft[i][0] *= window;
    fft[i][1] *= window;
  }

  fftwf_plan backward = fftwf_plan_dft_c2r_1d(in.size(), fft, out.data(), FFTW_ESTIMATE);
  fftwf_execute(backward);
  fftwf_free(fft);

  for (size_t i = 0; i < out.size(); i++)
  {
    out[i] /= out.size();
    out[i] += average;
  }

  return out;
}

std::vector<float> moving_average(const std::vector<float> &vector, const size_t period = 10){
  std::vector<float> output(vector.size() - period);
  for (size_t i = 0; i < vector.size() - period; i++){
    output[i] = std::accumulate(std::begin(vector) + i, std::begin(vector) + i + period, 0)/period;
  }
  return output;
}

size_t find_SOFAR_channel(const std::vector<float> &speed_of_sound)
{
  std::vector<std::vector<float>::const_iterator> minima;
  std::vector<std::vector<float>::const_iterator> maxima;
  std::vector<std::vector<float>::const_iterator> stationary;

  auto end = std::end(speed_of_sound);
  auto begin = 1 + std::begin(speed_of_sound);

  float minimum = MAXFLOAT;
  std::vector<float_t>::const_iterator min_pos = std::end(speed_of_sound);

  for (; begin != end - 1; begin++)
  {
    if (*begin < *(begin - 1) && *begin < *(begin + 1))
    {
      for (auto backtracker = begin; backtracker != std::begin(speed_of_sound); --backtracker)
      {
        if (*backtracker > *begin + 5)
        {
          for (auto tracker = begin; tracker != end; tracker++)
          {
            if (*tracker > *begin + 5)
            {
              minima.push_back(begin);
              if (minimum > *begin)
              {
                min_pos = begin;
                minimum = *begin;
              }
              break;
            }
          }

          break;
        }
      }
    }
    else if (*begin > *(begin - 1) && *begin > *(begin + 1))
    {
      maxima.push_back(begin);
    }
    else if (std::abs(*begin - *(begin - 1)) < 1e-4)
    {
      stationary.push_back(begin);
    }
  }

  return std::distance(std::begin(speed_of_sound), min_pos);
}

typedef std::vector<float> argodata_t;
typedef std::vector<argodata_t> argotable_t;

struct [[nodiscard]] argodata_struct {
  argotable_t depths;
  argotable_t speeds_of_sound;
  argodata_t lats;
  argodata_t longs;
  argodata_t dates;
};

const argodata_struct read_file_data(fs::path filepath){
  using namespace netCDF;
  NcFile datafile(filepath.string(), NcFile::read);

  static const size_t N_PROF = datafile.getDim("N_PROF").getSize();
  static const size_t N_LEVELS = datafile.getDim("N_LEVELS").getSize();

  argodata_struct data {
    argotable_t(),
    argotable_t(),
    argodata_t(N_PROF),
    argodata_t(N_PROF),
    argodata_t(N_PROF)
  };

  data.depths.reserve(N_PROF);
  data.speeds_of_sound.reserve(N_PROF);

  NcVar latsVar, longsVar, dateVar;
  latsVar = datafile.getVar("LATITUDE");
  longsVar = datafile.getVar("LONGITUDE");
  dateVar = datafile.getVar("JULD");

  NcVar tempVar, depthVar, salinityVar;
  tempVar = datafile.getVar("TEMP");
  depthVar = datafile.getVar("DEPH_CORRECTED");
  salinityVar = datafile.getVar("PSAL_CORRECTED");

  if (tempVar.isNull() || depthVar.isNull() || salinityVar.isNull())
    exit(1);

  latsVar.getVar(data.lats.data());
  longsVar.getVar(data.longs.data());
  dateVar.getVar(data.dates.data());


  for (size_t j = 0; j < N_PROF; j++)
  {
    auto tempIn = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
    auto salIn = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
    auto depthIn = std::vector<float_t>(N_LEVELS, std::nan("nodata"));
    auto speed_of_sound_vec =
        std::vector<float_t>(N_LEVELS, std::nan("nodata"));

    std::vector<size_t> start{j, 0};
    std::vector<size_t> count{1, N_LEVELS};

    salinityVar.getVar(start, count, salIn.data());
    tempVar.getVar(start, count, tempIn.data());
    depthVar.getVar(start, count, depthIn.data());

    for (size_t i = 0; i < N_LEVELS; i++)
    {
      if (tempIn[i] != 99999 && salIn[i] != 99999 && depthIn[i] != 99999)
      {
        speed_of_sound_vec[i] = speed_of_sound::speed_of_sound(
            speed_of_sound::pressure_at_depth(depthIn[i], data.lats[j]), tempIn[i],
            salIn[i]);
      }
      else
      {
        tempIn[i] = std::nan("nodata");
        salIn[i] = std::nan("nodata");
        depthIn[i] = std::nan("nodata");
      }
    }
    {

      auto it = std::begin(speed_of_sound_vec);
      for (; it != std::end(speed_of_sound_vec); ++it)
      {
        if (std::isnan(*it))
          break;
      }
      auto size = std::abs(std::distance(std::begin(speed_of_sound_vec), it));
      speed_of_sound_vec.resize(size);
      depthIn.resize(size);

      data.depths.push_back(depthIn);
      data.speeds_of_sound.push_back(speed_of_sound_vec);
    }
  }
  return std::move(data);
}

int main(int argc, char *argv[])
{

  std::cout << "ARGO (Met Office Hadley Centre) Data" << std::endl;

  if (argc == 1)
  {
    std::cout << "Please pass a filename to the programme" << std::endl;
    return 1;
  }

  fs::path argument (argv[1]);

  auto data = read_file_data(argument);
  

  std::cout << "[INFO]: FINISHED PROCESSING\n"
            << std::endl;

  std::fstream f("channels.txt", std::fstream::out);

  for (size_t i = 0; i < 10; i++)
  {
    auto sos_vect = data.speeds_of_sound[i];
    auto depth_vect = data.depths[i];
    auto index = find_SOFAR_channel(sos_vect);
    // f << "At " << lats[i] << ":" << longs[i] << " at " << dates[i];
    f  << ", the SOFAR minimum is " << depth_vect[index] << "m under the surface"
      << std::endl;

    /*
    std::vector<float_t> difference(N_LEVELS, 0);
    for (size_t j = 1; j < N_LEVELS; j++)
    {
            if (depth_vect[j] - depth_vect[j - 1] != 0 && !isnan(depth_vect[j])
    && !isnan(depth_vect[j-1])) difference[j] = (sos_vect[j] - sos_vect[j - 1])
    / (depth_vect[j] - depth_vect[j - 1]);
    }*/

    auto vect = low_pass_filter(sos_vect);
    depth_vect.resize(vect.size());

    Gnuplot gp;
    gp << "set multiplot layout 2,1\n"
       << "set key off\n";
    grapher::plot_lines(gp, depth_vect, vect, sos_vect);
    grapher::plot_lines(gp, moving_average(depth_vect), moving_average(sos_vect));
  }
  f.close();
  return 0;
}
