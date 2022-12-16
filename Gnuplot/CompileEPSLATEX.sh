#!/bin/bash
echo 'This Script compiles the epslatex in the standalone mode'
echo 'Syntax: $PATH/EPSLATEX2PNG.sh ["Gnuplot Makefile"] ["tex output"]'
echo '# SCRIPT WAS CALLED BY: #'
echo "${0} $@"
echo "Getting Standard Gnuplot output"
#set terminal  epslatex standalone color colortex 10
#Get output 
echo "Gnuplot Makefile: $1"
#echo "Latex output: $2"
latexout=`grep 'PLOTNAME="' $1 | tail -n 1 |  sed 's/PLOTNAME=//g' | sed 's/"//g'`

psfile=`echo "${latexout}.ps"`
dvifile=`echo "${latexout}.dvi"`
pngfile=`echo "${latexout}.png"`
#psfile=`echo "${2%.tex}.ps"`
#dvifile=`echo "${2%.tex}.dvi"`
#pngfile=`echo "${2%.tex}.png"`
echo "DVI output: $dvifile"
echo "Postscript output: $psfile"
echo "PNG output: $pngfile"

gnuplot "$1" 
latex "${latexout}.tex"
dvips -o $psfile $dvifile 
ps2pdf $psfile 
convert $psfile $pngfile 
