unset pm3d 
unset key
unset border
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


reset
set xrange [XMIN:XMAX]
set yrange [YMIN:YMAX]
set cbrange [CBMIN:CBMAX]   
print it, CBMIN, CBMAX, (CBMIN*abs(it)/1000.0), (CBMAX*abs(it)/1000.0)
CTINCREMENT=(CBMAX-CBMIN)*abs(it)/10000.0
CBMINSCALE=(CBMIN*abs(it)/1000.0)
CBMAXSCALE=(CBMAX*abs(it)/1000.0) 
print it,  CBMINSCALE+CTINCREMENT, CTINCREMENT, CBMAXSCALE-CTINCREMENT
#print CBMINSCALE+CTINCREMENT, (abs(it)/10000.0), CBMAXSCALE-CTINCREMENT


set isosample 250, 250
set table 'temp.dat'
set contour base
set cntr levels incremental CBMINSCALE+CTINCREMENT, CTINCREMENT, CBMAXSCALE-CTINCREMENT
set cntr bspline
unset surface
splot INPUT_FILE u 2:1:($3*Scaling_factor_Z) w line  lc palette
unset table

reset
unset border
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
set xrange [XMIN:XMAX]
set yrange [YMIN:YMAX]

set size ratio SIZE_RATIO
set lmargin at screen LMARGIN
set rmargin at screen RMARGIN
set tmargin at screen TMARGIN
set bmargin at screen BMARGIN
XPosition_Label= LMARGIN+ (0.2*(RMARGIN- LMARGIN))
YPosition_Label= TMARGIN- (0.2*(TMARGIN- BMARGIN))
LABEL=sprintf("%d mV",it);
set label LABEL at screen XPosition_Label, screen YPosition_Label  #font 'Helvetica,8'
plot 'temp.dat' w l lt -1 lw 0.5






