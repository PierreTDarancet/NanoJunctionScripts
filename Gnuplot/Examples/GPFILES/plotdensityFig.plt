unset pm3d
unset key
set colorbox
set surface
set border
unset title
unset contour
unset xlabel
unset ylabel
unset cblabel
unset xtics
unset ytics
unset cbtics

if (SETTITLE)   set title sprintf("%s Changes in Density", GENERIC_NAME)
set size ratio SIZE_RATIO
set lmargin at screen LMARGIN
set rmargin at screen RMARGIN
set tmargin at screen TMARGIN
set bmargin at screen BMARGIN
set xrange [XMIN:XMAX]

set yrange [YMIN+0.0001:YMAX]

set cbrange [CBMIN:CBMAX]   


if (SETXLABEL)   set xlabel XLABEL 
set xtics scale 2 XTICS 
#if (SETXLABEL) set xtics  XTICS  

if (SETYLABEL)   set ylabel YLABEL offset 0,7
set ytics scale 2 YTICS 

if (SETCBLABEL)  set cblabel CBLABEL offset 0,7
set cbtics scale 2 CBTICS
set colorbox

set palette defined (CBMIN "blue", 0 "white", CBMAX "red")
set palette defined (CBMIN "cyan", (CBMIN/2.) "blue", 0 "white", (CBMAX/2.0) "red", CBMAX "orange")
LABELPOSH=LMARGIN+0.06
LABELPOSV=TMARGIN-0.04
LABELNAME=""
LABELNAME=sprintf("%d mV", abs(it))


plot INPUT_FILE u 2:1:($3*Scaling_factor_Z)  lc palette

