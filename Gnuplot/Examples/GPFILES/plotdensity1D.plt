unset pm3d

set border
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

if (SETTITLE)   set title sprintf("%s Changes in Density", GENERIC_NAME)
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
set xtics  XTICS 
#if (SETXLABEL)   set xlabel at XOFFSET
if (SETYLABEL)   set ylabel YLABEL offset 0,7
set ytics  YTICS
#if (SETYLABEL)   set ylabel at YOFFSET  
if (SETCBLABEL)  set cblabel CBLABEL #CBOFFSET
if (SETCBLABEL)  set cbtics  CBTICS
#if (SETCBLABEL)  set cblabel at 

#print(Scaling_factor_Y)
#print(INPUT_FILE)
XPosition_Label= LMARGIN+ (0.2*(RMARGIN- LMARGIN))
YPosition_Label= TMARGIN- (0.2*(TMARGIN- BMARGIN))
LABEL=sprintf("%d mV",it);
set label LABEL at screen XPosition_Label, screen YPosition_Label  #font 'Helvetica,8'
plot INPUT_FILE u 1:($2*Scaling_factor_Y) with lines lw 3


