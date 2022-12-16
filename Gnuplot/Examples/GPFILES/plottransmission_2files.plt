unset pm3d
unset key
unset colorbox
unset surface
unset title 
unset contour
unset xlabel
unset ylabel
unset y2label

unset cblabel
unset xtics
unset ytics

unset x2tics
unset y2tics
unset cbtics

if (SETTITLE)   set title sprintf("%s: %d mV", GENERIC_NAME,it)
set size ratio SIZE_RATIO
set lmargin at screen LMARGIN
set rmargin at screen RMARGIN
set tmargin at screen TMARGIN
set bmargin at screen BMARGIN
set xrange [XMIN:XMAX]
XTICS=(XMAX-XMIN)/8
set yrange [(YMIN+0.00001):1]
set y2range [(YMIN+0.00001):1]

print SETYLABEL
if (SETXLABEL)   set xlabel XLABEL 
#set xtics  out
set xtics  XTICS font "Helvetica,0.1"
if (SETXLABEL) set xtics  XTICS   font 'Helvetica,10'
if (SETYLABEL)  set y2label YLABEL
set log y2
set y2tics  0, 100

LABELPOSH=LMARGIN
LABELPOSV=TMARGIN-0.03
LABELNAME=""
LABELNAME=sprintf("%d mV", it)

set label LABELNAME at screen LABELPOSH, screen LABELPOSV font 'Helvetica,10'

plot INPUT_FILE1 using ($1-FermiEnergy):2  axis x1y2 with filledcurve y2=-10 lt 5 lw 2,  INPUT_FILE1 u ($1-FermiEnergy):2   axis x1y2 with lines  lt -1 lw 1, INPUT_FILE2 using ($1-FermiEnergy):2  axis x1y2 with filledcurve y2=-10 lt 5 lw 2,  INPUT_FILE2  u ($1-FermiEnergy):2   axis x1y2 with lines  lt -1 lw 1, INPUT_FILE3  u ($1-FermiEnergy):2  axis x1y2 with impulses lt 1 lw 2
