set terminal png font '/usr/share/fonts/newpx/TeXGyrePagellaX-Regular.otf,12'
set datafile separator ','

system("mkdir -p output/figures")
speeds_of_sound = system('ls output/speed_of_sound/*.csv')
temperatures = system('ls output/temp/*.csv')

f(x) = a*x + b
set key autotitle columnhead
do for [file in speeds_of_sound] {
    filename = "./output/figures/sos-" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set yrange [2000:0]
    plot file u 1:2;
}

do for [file in temperatures] {
    filename = "./output/figures/temp-" . system("basename " . file . " .csv") . ".png"
    print filename
    set output filename
    set yrange [2000:0]
    plot file u 1:2;
}