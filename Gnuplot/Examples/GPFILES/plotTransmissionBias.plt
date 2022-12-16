reset
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


set style line 1 lw 2 lt 1 lc rgb "orange"
set style line 2 lw 2 lt 2 lc rgb "red"
set style line 3 lw 3 lt 1 lc rgb "black"
set style line 4 lw 4 lt 3 lc rgb "blue"
set style line 5 lw 4 lt 2 lc rgb "cyan"

set size ratio SIZE_RATIO
set lmargin at screen LMARGIN
set rmargin at screen RMARGIN
set tmargin at screen TMARGIN
set bmargin at screen BMARGIN
set xrange [XMIN:XMAX]
set yrange [YMIN:YMAX]


#print(XLABEL);
if (SETXLABEL)   set xlabel XLABEL font 'Helvetica,20' offset 0,0.8
set xtics  XTICS font 'Helvetica,0.01'
set mxtics  5 
if (SETXLABEL)  set xtics  XTICS font 'Helvetica,20' offset 0,0.2

#if (SETXLABEL)   set xlabel at XOFFSET
if (SETYLABEL)   set ylabel YLABEL font 'Helvetica,20' offset 0,6
set ytics  YTICS font 'Helvetica,20' offset 0.2,0


unset key 

if (SETCBLABEL) set label GENERIC_NAME at screen  0.8, screen 0.95 font 'Helvetica,20'
if (SETCBLABEL) set key at screen HKEY, screen VKEY


BIAS1="-1000"
BIAS2="-500"
BIAS3="0"
BIAS4="500"
BIAS5="1000"
#LABEL1=sprintf("%s mV",BIAS1);
#LABEL2=sprintf("%s mV",BIAS2);
#LABEL3=sprintf("%s mV",BIAS3);
#LABEL4=sprintf("%s mV",BIAS4);
#LABEL5=sprintf("%s mV",BIAS5);
LABEL1="-1V"
LABEL2="-0.5V"
LABEL3="0V"
LABEL4="0.5V"
LABEL5="1V"



INPUT_FILE_1_A=sprintf("%s.%s",BIAS1, INPUT_FILE);
INPUT_FILE_2_A=sprintf("%s.%s",BIAS2, INPUT_FILE);
INPUT_FILE_3_A=sprintf("%s.%s",BIAS3, INPUT_FILE);
INPUT_FILE_4_A=sprintf("%s.%s",BIAS4, INPUT_FILE);
INPUT_FILE_5_A=sprintf("%s.%s",BIAS5, INPUT_FILE);



set log y
set format y "10^{%L}"
set yrange [YMIN:YMAX]
set ytics 10 font 'Helvetica,20'
set mytics 10 
print(YMIN)

set label LABEL at screen  LMARGIN+HLABEL, TMARGIN+VLABEL font 'Helvetica,20'
plot INPUT_FILE_1_A using ($1-FermiEnergy):2  w l ls 1 title LABEL1,\
     INPUT_FILE_2_A using ($1-FermiEnergy):2  w l ls 2 title LABEL2,\
     INPUT_FILE_3_A using ($1-FermiEnergy):2  w l ls 3 title LABEL3,\
     INPUT_FILE_4_A using ($1-FermiEnergy):2  w l ls 4 title LABEL4,\
     INPUT_FILE_5_A using ($1-FermiEnergy):2  w l ls 5 title LABEL5




