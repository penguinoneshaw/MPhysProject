set terminal png font '/usr/share/fonts/newpx/TeXGyrePagellaX-Regular.otf,12'
set datafile separator ','

system("mkdir -p output/figures")
speeds_of_sound = system('ls output/speed_of_sound/*.csv')
temperatures = system('ls output/temp/*.csv')
power_spectra = system('ls output/power_spectra/*.csv')


f(x) = a*x + b
set key autotitle columnhead
set key off
do for [file in temperatures] {
    filename = "./output/figures/temp-" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set xlabel "Time (days since 1950-01-01)"
    set ylabel "Temperature at predicted minimum"
    fit f(x) file u 1:2 via a,b;

    plot file u 1:2, f(x);
}

do for [file in power_spectra] {
    filename = "./output/figures/ps-" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set xlabel "Frequency"
    set ylabel "Normalised Power Spectrum"

    plot file u ($3/3000) w lines;
}

do for [file in speeds_of_sound] {
    filename = "./output/figures/sos-" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set xlabel "Time (days since 1950-01-01)"
    set ylabel "Depth below sea level (m)"
    fit f(x) file u 1:2 via a,b;
    set yrange [2000:0]
    plot file u 1:2, f(x);
}
