/** @file grapher.h
 * Defines functions for graphing vectors of samples from a PDF.
 */

#if !defined(H_GRAPHER)
#define H_GRAPHER

#include "../gnuplot-iostream/gnuplot-iostream.h"
#include <vector>
#include <iostream>
#include <valarray>

namespace grapher
{
  /** @namespace grapher
   * Clusters together the utility functions for plotting.
   */

  template <typename T>
  void histogram_distributions(std::vector<T> points, std::string title = "Probability Distribution")
  {
    /**
     * Uses Gnuplot to generate a histogram for a given set of points.
     * Uses the Freedman-Diaconis rule to set the bin widths and generates the histogram on the Gnuplot side.
     *
     * @param points A vector of the points to be histogrammed.
     * @param title A string which will become the plot title.
     */
    std::sort(std::begin(points), std::end(points));

    T IQR = points[std::floor(points.size() * 3 / 4)] - points[std::floor(points.size() * 1 / 4)];

    T binwidth = 2 * IQR * std::pow(points.size(), -1 / 3.);
    T xmin = points[0];
    T xmax = points[points.size() - 1];

    Gnuplot gp;
    gp << "set xrange [" << xmin << ":" << xmax << "]\n";
    gp << "set style fill solid 0.5\n";
    gp << "binwidth = " << binwidth << "\nbinstart = " << xmin << "\nset boxwidth binwidth\n";
    gp << "plot '-' using (binwidth*(floor(($1-binstart)/binwidth)+0.5)+binstart):(1.0) smooth freq w boxes title '" << title << "'\n";

    gp.send1d(points);
  }

  template <typename T, typename F>
  void plot_points(T const &xs, F const &ys)
  {
    /**
     * Scatter plots a pair of vectors of samples.
     *
     * @param xs Vector of x values
     * @param ys Vector of y values
     */
    Gnuplot gp;
    // gp << "set xrange [0:10]\nset yrange[0:" << 2*M_PI << "]\n";
    gp << "plot '-' with points pt 4\n";

    gp.send1d(std::make_tuple(xs, ys));
  }

  template <typename F, typename T>
  int plot_line(Gnuplot &gp, F const &ys, const T &xs)
  {
    /**
     * Scatter plots a pair of vectors of samples.
     *
     * @param gp Reference to the gnuplot connection object
     * @param xs Vector of x values
     * @param ys Vector of y values
     */
    gp << gp.file1d(std::make_tuple(xs, ys)) << "with lines ls 2, \\\n";
   
    return 0;
  }

  template <typename F, typename ... Ts>
  void plot_lines(Gnuplot &gp, F const &ys, Ts ... xss){
    gp << "plot ";
    [](...){}(plot_line(gp, ys, xss)...);
    gp << std::endl;
  }
}

#endif // H_GRAPHER
