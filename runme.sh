#!/usr/bin/bash

clean()
{
if test -f ./StVector.txt;then
    rm ./StVector.txt
fi
if test -e ./run;then
    rm ./run
fi
if test -f ./Orbit6elem.txt;then 
    rm ./Orbit6elem.txt
fi
if test -f ./a.svg;then
    rm ./*.svg
fi
}
draw()
{
gnuplot << WTF
set terminal wxt size 1920,1080 enhanced font 'Verdana,10' persist
unset key
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set xyplane at 0
set hidden3d
set xrange [-8e6:8e6]
set yrange [-8e6:8e6]
set zrange [-8e6:8e6]
set view 45,45
set view equal xyz
set title "Orbit"
R=6371004
f(x,y)=(R**2-x**2-y**2)**0.5
g(x,y)=-(R**2-x**2-y**2)**0.5
set isosamples 30
set grid
splot f(x,y),\
      g(x,y),\
      'StVector.txt' using 5:6:7 with lines linecolor rgb "black" linewidth 1.75
WTF
gnuplot << WTF
    unset key
    set terminal svg size 1280,720 enhanced font 'Verdana,10'
    set output "a.svg"  
    set object 1 rectangle from screen 0,0 to screen 1,1 fc rgb "white10" behind
    set xlabel "time(s)"
    set ylabel "a(m)"
    set title "a"
    set grid
    plot\
    'Orbit6elem.txt' using 1:2 with lines\
    linecolor 3 linewidth 1.75
WTF
gnuplot << WTF
    unset key
    set terminal svg size 1280,720 enhanced font 'Verdana,10'
    set output "e.svg"  
    set object 1 rectangle from screen 0,0 to screen 1,1 fc rgb "white10" behind
    set xlabel "time(s)"
    set ylabel "e"
    set title "e"
    set grid
    plot\
    'Orbit6elem.txt' using 1:3 with lines\
    linecolor 3 linewidth 1.75 axis x1y1
WTF
gnuplot << WTF
    unset key
    set terminal svg size 1280,720 enhanced font 'Verdana,10'
    set output "i.svg"  
    set object 1 rectangle from screen 0,0 to screen 1,1 fc rgb "white10" behind
    set xlabel "time(s)"
    set ylabel "i(rad)"
    set title "i"
    set grid
    plot\
    'Orbit6elem.txt' using 1:4 with lines\
    linecolor 3 linewidth 1.75 axis x1y1
WTF
gnuplot << WTF
    unset key
    set terminal svg size 1280,720 enhanced font 'Verdana,10'
    set output "O.svg"  
    set object 1 rectangle from screen 0,0 to screen 1,1 fc rgb "white10" behind
    set xlabel "time(s)"
    set ylabel "O(rad)"
    set title "O"
    set grid
    plot\
    'Orbit6elem.txt' using 1:5 with lines\
    linecolor 3 linewidth 1.75 axis x1y1
WTF
gnuplot << WTF
    unset key
    set terminal svg size 1280,720 enhanced font 'Verdana,10'
    set output "w.svg"  
    set object 1 rectangle from screen 0,0 to screen 1,1 fc rgb "white10" behind
    set xlabel "time(s)"
    set ylabel "w(rad)"
    set title "w"
    set grid
    plot\
    'Orbit6elem.txt' using 1:6 with lines\
    linecolor 3 linewidth 1.75 axis x1y1
WTF
gnuplot << WTF
    unset key
    set terminal svg size 1280,720 enhanced font 'Verdana,10'
    set output "f.svg"  
    set object 1 rectangle from screen 0,0 to screen 1,1 fc rgb "white10" behind
    set xlabel "time(s)"
    set ylabel "f(rad)"
    set title "f"
    set grid
    plot\
    'Orbit6elem.txt' using 1:7 with lines\
    linecolor 3 linewidth 1.75 axis x1y1
WTF
}
generate()
{
gnuplot << WTF
set terminal svg size 1920,1080 enhanced font 'Verdana,10' background "white10"
set output "orbit.svg"
unset key
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set xyplane at 0
set hidden3d
set xrange [-7e6:7e6]
set yrange [-7e6:7e6]
set zrange [-7e6:7e6]
set view $1,$2
set view equal xyz
set title "Orbit"
R=6371004
f(x,y)=(R**2-x**2-y**2)**0.5
g(x,y)=-(R**2-x**2-y**2)**0.5
set isosamples 30
set grid
splot f(x,y),\
      g(x,y),\
      'data' using 5:6:7 with lines linecolor rgb "black" linewidth 1.75
WTF
}
main()
{
if test -f ./sat.c;then
    if test ! -f ./run;then
	gcc -lm -lgsl -lcblas sat.c -o run
    fi
    ./run
fi
draw
}
if test $# != 0;then
    if test "$1" == "clean";then
        echo cleaning...
        clean
        echo "clean done!"
    elif test "$1" == "gene";then
        echo "generating svg file..."
        generate "$2" "$3"
        echo "generating svg file done!"
    elif "$1" == "draw";then
        draw
    fi
else
    main
fi

