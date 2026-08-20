#include <bits/chrono.h>
#include <time.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <concepts>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <netcdf>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "fitting/projectfit.hpp"
#include "ncChar.h"
#include "ncException.h"
#include "oQTM/oQTM.hpp"
#include "speed_of_sound.hpp"

namespace fs = std::filesystem;

namespace project {

template <typename T, typename K, std::floating_point V>
std::tuple<std::vector<std::tuple<T, T, K, V>>,
           std::vector<std::tuple<T, T, K, V>>,
           std::vector<std::tuple<T, T, K, V>>>
read_to_tree(const fs::path& filepath,
             speed_of_sound::SpeedOfSoundAlgorithm method =
                 speed_of_sound::ALGORITHM_UNESCO,
             const fit::FitFunction fitfunction = fit::FIT_QUADRATIC) {
    using namespace netCDF;

    NcFile datafile(filepath.string(), NcFile::read);
    std::cout << filepath.string() << std::endl;
    static const std::size_t N_PROF = datafile.getDim("N_PROF").getSize();
    static const std::size_t N_LEVELS = datafile.getDim("N_LEVELS").getSize();

    std::vector<std::tuple<T, T, K, V>> results;
    std::vector<std::tuple<T, T, K, V>> temp_results;
    std::vector<std::tuple<T, T, K, V>> error_results;

    // results.reserve(N_PROF);

    NcVar latsVar, longsVar, dateVar, platformVar, tempVar, depthVar, salinityVar,
        profileQualityVar, levelsQualityVar, positionQualityVar, tempQualityVar;
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

    std::vector<double> lats(latsVar.getDim(0).getSize(), 0),
        lons(longsVar.getDim(0).getSize(), 0);
    std::vector<double> dates(dateVar.getDim(0).getSize(), 0);

    if (tempVar.isNull() || depthVar.isNull() || salinityVar.isNull()) {
        throw std::runtime_error("File format not as expected.");
    }

    latsVar.getVar(lats.data());
    longsVar.getVar(lons.data());
    dateVar.getVar(dates.data());

    for (std::size_t i = 0; i < N_PROF; i++) {
        // Iterate over profiles
        try {
            uint32_t profile_qc;

            char position_qc = 0;
            char potm_qc = 0;

            profileQualityVar.getVar(std::vector<size_t>{i}, &profile_qc);
            positionQualityVar.getVar(std::vector<size_t>{i}, &position_qc);
            tempQualityVar.getVar(std::vector<size_t>{i}, &potm_qc);

            if ((profile_qc & 0b11) != 0 || position_qc == 0b11 || potm_qc == 0b11) {
                // Skip if the profile fails a quality check — see the datasource for
                // details
                continue;
            }
        } catch (const netCDF::exceptions::NcInvalidCoords& e) {
            break;
        }

        T lat = lats[i];
        T lon = lons[i];
        K date = dates[i];

        // reserve space in all data for a profile
        std::vector<float_t> tempIn(tempVar.getDim(0).getSize(),
                                    std::nan("nodata"));
        std::vector<float_t> salIn(salinityVar.getDim(0).getSize(),
                                   std::nan("nodata"));
        std::vector<float_t> depthIn(depthVar.getDim(0).getSize(),
                                     std::nan("nodata"));
        std::vector<uint32_t> level_qc_bits(levelsQualityVar.getDim(0).getSize(),
                                            std::nan("nodata"));

        std::vector<double> actual_depths, actual_temperatures;
        std::vector<double> speed_of_sound_vec;

        // Read in data
        actual_depths.reserve(N_LEVELS);
        speed_of_sound_vec.reserve(N_LEVELS);
        actual_temperatures.reserve(N_LEVELS);
        std::vector<std::size_t> start{i, 0};
        std::vector<std::size_t> temp_count{1, tempVar.getDim(1).getSize()};
        std::vector<std::size_t> salinity_count{1, salinityVar.getDim(1).getSize()};
        std::vector<std::size_t> depth_count{1, depthVar.getDim(1).getSize()};
        std::vector<std::size_t> levels_quality_count{
            1, levelsQualityVar.getDim(1).getSize()};

        try {
            levelsQualityVar.getVar(start, levels_quality_count,
                                    level_qc_bits.data());
            salinityVar.getVar(start, salinity_count, salIn.data());
            tempVar.getVar(start, temp_count, tempIn.data());
            depthVar.getVar(start, depth_count, depthIn.data());

        } catch (const netCDF::exceptions::NcInvalidCoords& e) {
            std::cerr << filepath << ": " << e.what() << std::endl;
            break;
        }

        // Process to a speed of sound profile
        for (std::size_t j = 0; j < N_LEVELS; j++) {
            if ((level_qc_bits[j]) != 0 || std::isnan(level_qc_bits[j]))
                continue;
            T s{};
            // This is probably templatable
            switch (method) {
                case speed_of_sound::ALGORITHM_LEROY:
                    s = speed_of_sound::leroy_et_al((double)depthIn[j], (double)tempIn[j],
                                                    (double)salIn[j], lat);
                    break;

                case speed_of_sound::ALGORITHM_UNESCO:
                    s = speed_of_sound::speed_of_sound(
                        speed_of_sound::pressure_at_depth((double)depthIn[j], lat),
                        (double)tempIn[j], (double)salIn[j]);
                    break;
            }

            speed_of_sound_vec.push_back(s);
            actual_depths.push_back((double)depthIn[j]);
            actual_temperatures.push_back(tempIn[j]);
        }

        // Validate that there's enough data to continue
        if (actual_depths.size() < 10) {
            continue;
        };

        try {
            V xmin, errmin;
            switch (fitfunction) {
                case fit::FIT_IDEAL:
                    std::tie(xmin, errmin) = fit::find_sofar_channel<V, fit::FIT_IDEAL>(
                        speed_of_sound_vec, actual_depths, 3);
                    break;

                default:
                    std::tie(xmin, errmin) = fit::find_sofar_channel<V, fit::FIT_QUADRATIC>(
                        speed_of_sound_vec, actual_depths, 3);
                    break;
            }
            if (xmin > actual_depths.back()) {
                continue;
            }
            size_t distance =
                std::distance(actual_depths.begin(),
                              std::find_if(actual_depths.begin(), actual_depths.end(),
                                           [=](auto a) { return a > xmin; }));
            auto minimum_temperature =
                actual_temperatures[distance >= (size_t)actual_temperatures.size()
                                        ? actual_temperatures.size() - 1
                                        : distance];

            auto surface_temperature = actual_temperatures[0];
            auto temperature_difference = minimum_temperature - surface_temperature;

            // auto tavg = std::accumulate(actual_temperatures.begin(),
            //                            actual_temperatures.end(), 0) /
            //            ((V)actual_temperatures.size());

            const V result{xmin};
            results.push_back(std::tie(lon, lat, date, result));

            const V error{errmin};
            error_results.push_back(std::tie(lon, lat, date, error));

            if (temperature_difference < -40 || temperature_difference > 40.0)
                continue;
            const V result_temp{temperature_difference};

            temp_results.push_back(std::tie(lon, lat, date, result_temp));
        } catch (const std::runtime_error& e) {
            continue;
        }
    }

    return std::tie(results, temp_results, error_results);
}
}  // namespace project

using namespace std::chrono_literals;
constexpr auto BASE_DATE =
    std::chrono::sys_days{1950y / std::chrono::January / 1d};

int main(int argc, char* argv[]) {
    using namespace project;
    typedef double value_t;
    typedef oQTM_Mesh<double, uint64_t, value_t, 6> mesh_t;

    mesh_t globemesh;
    mesh_t tempmesh;
    oQTM_Mesh<double, uint64_t, value_t, 6, true> errormesh;
    std::filesystem::path data_dir{"./data"};
    if (argc >= 2) {
        data_dir = std::filesystem::path{argv[1]};
    }

    speed_of_sound::SpeedOfSoundAlgorithm algorithm =
        speed_of_sound::ALGORITHM_UNESCO;

    fit::FitFunction fitfunction = fit::FIT_QUADRATIC;
    try {
        if (argc >= 3) {
            std::string options(argv[2]);
            if (options.find('l') != std::string::npos)
                algorithm = speed_of_sound::ALGORITHM_LEROY;
            if (options.find('i') != std::string::npos)
                fitfunction = fit::FIT_IDEAL;
        }
    } catch (const std::logic_error& e) {
        std::cerr << e.what() << '\n';
    }

    uint64_t GRANULARITY = 10;

    if (argc >= 4) {
        GRANULARITY = std::stoi(argv[3]);
    }

    auto dirs = fs::directory_iterator(data_dir);

    // #pragma omp parallel
    for (const auto& path : dirs) {
        if (fs::is_regular_file(path) && path.path().extension() == ".nc") {
            auto [points, temp_points, error_points] =
                read_to_tree<double, uint64_t, value_t>(path, algorithm, fitfunction);

            globemesh.insert(points);
            tempmesh.insert(temp_points);
            errormesh.insert(error_points);
        }
    }

    std::cout << "[INFO]: Interpolation complete" << std::endl;

    std::vector<std::pair<double, double>> locations{
        {-40.671900, 30.000},     // North Atlantic
        {-14.442143, -30.000},    // South Atlantic
        {18.179810, 35.013669},   // Mediterranean
        {79.738750, -21.668352},  // Indian Ocean
        {114.008616, 15.425780},  // South China Sea
        {158.790808, 14.856203},  // Western Pacific
        {-105.585039, 9.403472},  // Eastern Pacific
        {3.343785, 56.3836}       // North Sea
    };

    const std::string OUTPUT_DIRECTORY =
        "output-" + std::to_string(algorithm) + "-" +
        std::to_string(fitfunction) + "-" +
        std::to_string(
            std::chrono::system_clock::now().time_since_epoch().count() /
            86400000);

    fs::create_directory(OUTPUT_DIRECTORY);
    fs::create_directory(OUTPUT_DIRECTORY + "/power_spectra");
    fs::create_directory(OUTPUT_DIRECTORY + "/speed_of_sound");
    fs::create_directory(OUTPUT_DIRECTORY + "/temperature");

    for (uint8_t j = 0; j < 8; j++) {
        auto loc = mesh_t::location_t{j};
        auto averages = globemesh.get_averaged_points(loc, 1, GRANULARITY);
        auto temperatures = tempmesh.get_averaged_points(loc, 1, GRANULARITY);
        auto errors = errormesh.get_averaged_points(loc, 1, GRANULARITY);
        {
            std::stringstream filename;
            filename << OUTPUT_DIRECTORY + "/speed_of_sound/";

            filename << (std::size_t)j << ".sos-results.csv";
            std::ofstream fileout(filename.str());
            fileout
                << "days since 1950-01-01,date,speed of sound minimum depth (m),err"
                << std::endl;

            for (auto i : averages) {
                auto t = BASE_DATE + std::chrono::days(i.first);
                fileout << i.first << ',' << t << "," << i.second << ','
                        << errors[i.first] << std::endl;
            }
            fileout.close();
        }

        {
            std::stringstream filename;
            filename << OUTPUT_DIRECTORY + "/temperature/";

            filename << (std::size_t)j << ".temp-results.csv";
            std::ofstream fileout(filename.str());
            fileout << "days since 1950-01-01,date,sea surface temperatures"
                    << std::endl;

            for (auto i : temperatures) {
                auto t = BASE_DATE + std::chrono::days(i.first);

                fileout << i.first << ',' << t << "," << i.second << std::endl;
            }

            fileout.close();
        }

        try {
            auto power_spectrum = fit::analyse_periodicity(averages, GRANULARITY);
            std::stringstream filename_ps;

            filename_ps << OUTPUT_DIRECTORY + "/power_spectra/";
            filename_ps << (std::size_t)loc[0] << ".ps-results.csv";
            std::ofstream fileout(filename_ps.str());
            std::vector<double> absolutes(power_spectrum.size());
            std::transform(power_spectrum.begin(), power_spectrum.end(),
                           absolutes.begin(),
                           [](auto a) { return std::abs(a.second); });

            auto max = std::max_element(absolutes.begin(), absolutes.end());
            fileout << "frequency,Real,Imag,Normalised Power,Phase,Power"
                    << std::endl;

            for (auto [f, i] : power_spectrum) {
                fileout << f << "," << i.real() << "," << i.imag() << ","
                        << std::abs(i) / *max << "," << std::arg(i) << ","
                        << std::abs(i) << std::endl;
            }
            fileout.close();
        } catch (const std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
            continue;
        }
    }

    std::ofstream file("location_dict.csv");
    for (const auto& location : locations) {
        auto loc = globemesh.location(location.first, location.second);
        std::stringstream location_string;

        for (auto i : loc)
            location_string << (int)i;

        file << location_string.str() << "," << location.second << ","
             << location.first << std::endl;

        std::vector<int> mesh_depths{3, 6};

        for (int k : mesh_depths) {
            auto a = globemesh.get_averaged_points(loc, k, GRANULARITY);
            auto t = tempmesh.get_averaged_points(loc, k, GRANULARITY);
            auto errors = errormesh.get_averaged_points(loc, k, GRANULARITY);

            {
                std::ofstream fileout(OUTPUT_DIRECTORY + "/speed_of_sound/" +
                                      location_string.str() + "." + std::to_string(k) +
                                      ".sos-results.csv");
                fileout
                    << "days since 1950-01-01,date,speed of sound minimum depth (m),err"
                    << std::endl;

                for (auto i : a) {
                    auto t = BASE_DATE + std::chrono::days(i.first);
                    fileout << i.first << ',' << t << "," << i.second << ","
                            << errors[i.first] << std::endl;
                }
                fileout.close();
            }

            {
                std::ofstream fileout(OUTPUT_DIRECTORY + "/temp/" +
                                      location_string.str() + "." + std::to_string(k) +
                                      ".temp-results.csv");

                fileout << "days since 1950-01-01,date,sea surface temperatures"
                        << std::endl;
                for (auto i : t) {
                    auto t = BASE_DATE + std::chrono::days(i.first);
                    fileout << i.first << ',' << t << "," << i.second << std::endl;
                }
                fileout.close();
            }

            try {
                auto power_spectrum = fit::analyse_periodicity(a, GRANULARITY);
                std::ofstream fileout(OUTPUT_DIRECTORY + "/power_spectra/" +
                                      location_string.str() + "." + std::to_string(k) +
                                      ".ps-results.csv");
                auto absolutes = power_spectrum | std::views::transform([](auto a) {
                                     return std::abs(a.second);
                                 });
                auto max = std::max_element(absolutes.begin(), absolutes.end());
                fileout << "frequency,Real,Imag,Normalised Power,Phase,Power"
                        << std::endl;

                // unsigned int period = a.end()->first - a.begin()->first;
                for (auto [f, i] : power_spectrum) {
                    fileout << f << "," << i.real() << "," << i.imag() << ","
                            << std::abs(i) / *max << "," << std::arg(i) << ","
                            << std::abs(i) << std::endl;
                }

                fileout.close();
            } catch (const std::runtime_error& e) {
                std::cerr << e.what() << std::endl;
                continue;
            }
        }
    }
    file.close();
}