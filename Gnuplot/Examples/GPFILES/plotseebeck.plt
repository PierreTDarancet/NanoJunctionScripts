unset pm3d
unset key
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

#if (SETYLABEL)   set ylabel at YOFFSET  
#if (SETCBLABEL)  set cblabel CBLABEL #CBOFFSET
#if (SETCBLABEL)  set cbtics  CBTICS
#if (SETCBLABEL)  set cblabel at 
set yrange [*:*]
if (SETYLABEL)   set ytics  #YTICS
#print(Scaling_factor_Y)
#print(INPUT_FILE)
plot seebeckcoefficient(x)  with filledcurve y1=0  lt 4 lw 2, seebeckcoefficient(x)  lt 1 lw 1

