set terminal png font '/usr/share/fonts/newpx/TeXGyrePagellaX-Regular.otf,12'
set datafile separator ','
set timefmt "%Y-%m-%d"
system("mkdir -p figures")

speeds_of_sound = system('ls speed_of_sound/*.csv')
temperatures = system('ls temp/*.csv')
power_spectra = system('ls power_spectra/*.csv')

f(x) = a*x + b
g(x) = f*x*x + g * x + h
set key autotitle columnhead
set key off

do for [file in power_spectra] {
    filename = "./figures/" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set xlabel "Frequency"
    set ylabel "Normalised Power Spectrum"

    plot file u 1:4 w lines;
}

do for [file in temperatures] {
    filename = "./figures/" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set xdata time
    set format x "%Y"
    set xlabel "Date"
    set ylabel "Mean Temperature Across Profiles (degree Celsius)"
    fit f(x) file u 2:3 via a,b;
    fit g(x) file u 2:3 via f,g,h;
    plot file u 2:3, f(x), g(x);
}

do for [file in speeds_of_sound] {
    filename = "./figures/" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set xdata time
    set format x "%Y"
    set xlabel "Date"
    set ylabel "Depth below sea level (m)"
    fit f(x) file u 2:3:4 zerror via a,b;
    fit g(x) file u 2:3:4 zerror via f,g,h;

    set yrange [2000:0]
    plot file u 2:3:4 w yerr, f(x),g(x);
}
