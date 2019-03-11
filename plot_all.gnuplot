set terminal png font '/usr/share/fonts/newpx/TeXGyrePagellaX-Regular.otf,12'
set datafile separator ','
set timefmt "%Y-%m-%d"
system("mkdir -p figures")

speeds_of_sound = system('ls speed_of_sound/*.csv')
temperatures = system('ls temp/*.csv')
power_spectra = system('ls power_spectra/*.csv')

f(x) = a*x + b
set key autotitle columnhead
set key off

do for [file in power_spectra] {
    filename = "./figures/ps-" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set xlabel "Frequency"
    set ylabel "Normalised Power Spectrum"

    plot file u 0:3 w lines;
}

do for [file in temperatures] {
    filename = "./figures/temp-" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set xdata time
    set format x "%Y-%m-%d"
    set xlabel "Time (days since 1950-01-01)"
    set ylabel "Mean Temperature Across Profile (degree Celsius)"
    fit f(x) file u 3:2 via a,b;

    plot file u 3:2, f(x);
}

do for [file in speeds_of_sound] {
    filename = "./figures/sos-" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set xdata time
    set format x "%Y-%m-%d"
    set xlabel "Time (days since 1950-01-01)"
    set ylabel "Depth below sea level (m)"
    fit f(x) file u 3:2 via a,b;
    set yrange [2000:0]
    plot file u 3:2, f(x);
}
