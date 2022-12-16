set pm3d map

unset colorbox
unset surface
unset title 
unset contour
unset xlabel
unset ylabel
unset cblabel
unset xtics
unset ytics
unset cbtics

if (SETTITLE)   set title sprintf("%s: %d mV", GENERIC_NAME,it)
set size ratio SIZE_RATIO
set lmargin at screen LMARGIN
set rmargin at screen RMARGIN
set tmargin at screen TMARGIN
set bmargin at screen BMARGIN
set xrange [XMIN:XMAX]

set yrange [YMIN:YMAX]

set cbrange [CBMIN:CBMAX]   

#print(XLABEL);
if (SETXLABEL)   set xlabel XLABEL #XOFFSET
if (SETXLABEL)   set xtics  XTICS 
#if (SETXLABEL)   set xlabel at XOFFSET
if (SETYLABEL)   set ylabel YLABEL #YOFFSET
if (SETYLABEL)   set ytics  YTICS
#if (SETYLABEL)   set ylabel at YOFFSET  
if (SETCBLABEL)  set cblabel CBLABEL #CBOFFSET
if (SETCBLABEL)  set cbtics  CBTICS
#if (SETCBLABEL)  set cblabel at 

#print(Scaling_factor_Y)
#print(INPUT_FILE)
#plot INPUT_FILE  u ($1+FermiEnergy):2 with filledcurve y1=0  lt 5 lw 2, INPUT_FILE u ($1+FermiEnergy):2  lt -1 lw 1
plot fitfunction(x) u 1:2 with points  2
