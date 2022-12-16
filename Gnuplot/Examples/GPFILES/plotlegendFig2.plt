#########################################
########       Color Box         ########
#########################################

set  border
set pm3d map
set size ratio HeightColorbox/WidthColorbox
LMARGIN = LMARGIN_CB
RMARGIN = RMARGIN_CB 
TMARGIN = TMARGIN_CB
BMARGIN = BMARGIN_CB
set lmargin at screen LMARGIN
set rmargin at screen RMARGIN
set tmargin at screen TMARGIN
set bmargin at screen BMARGIN
unset title
unset contour
unset colorbox
unset xlabel
unset ylabel
unset cblabel
unset xtics
unset ytics
unset cbtics


set pm3d map
g(y)=y
unset key
set xrange [-50:50]
set yrange [-50:50]
#set y2range [-50:50]
set zrange [-50:50]
set cbrange [-50:50]
set ylabel "Difference in Hartree Potential [V]"  font 'Helvetica,10'
set ytics ("-V/2" -50, " " -40 , " " -30 , " " -20 , " " -10 ,"0" 0, " " 10,  " " 20, " " 30, " " 40, "V/2" 50)

set isosamples 101, 101
set palette defined (-50 "blue", 0 "white", 50 "red")
splot g(y) 

#set y2label "Difference in Charge Density [e/]"  font 'Helvetica,10'
#set y2tics ("Check" -50, "0" 0,"Check" 50)

