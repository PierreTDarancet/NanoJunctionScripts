set pm3d map

unset key
set colorbox
set surface
unset title
unset contour
unset xlabel
unset ylabel
unset cblabel
unset xtics
unset ytics
unset cbtics

if (SETTITLE)   set title sprintf("%s: %d mV", GENERIC_NAME, it)
set size ratio SIZE_RATIO
set lmargin at screen LMARGIN
set rmargin at screen RMARGIN
set tmargin at screen TMARGIN
set bmargin at screen BMARGIN
set xrange [XMIN:XMAX]

set yrange [YMIN:YMAX]

set cbrange [CBMIN:CBMAX]   

if (SETXLABEL)   set xlabel XLABEL #XOFFSET
if (SETXLABEL)   set xtics  XTICS 
#if (SETXLABEL)   set xlabel at XOFFSET
if (SETYLABEL)   set ylabel YLABEL #YOFFSET
if (SETYLABEL)   set ytics  YTICS
#if (SETYLABEL)   set ylabel at YOFFSET  
if (SETCBLABEL)  set cblabel CBLABEL #CBOFFSET
if (SETCBLABEL)  set cbtics  CBTICS
#if (SETCBLABEL)  set cblabel at 

set palette defined (CBMIN "blue", 0 "white", CBMAX "red")
#print(Scaling_factor_Z)
#print(INPUT_FILE)
splot INPUT_FILE u 2:1:($3*Scaling_factor_Z)


