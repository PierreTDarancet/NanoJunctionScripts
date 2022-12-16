

~/SRC/ScriptsSiesta/Raman/genconstr.exe 4 28


Name="Si2Ch05" ; sourcedir=`pwd` ; for dir in $(ls Plus_* Minus_* ) ; do mkdir tmp_${dir} ; cd ${sourcedir}/tmp_${dir} ; pwd ; cp   ${sourcedir}/*.psf  ${sourcedir}/tmp_${dir}/ ;  cp   ${sourcedir}/pbsscript   ${sourcedir}/tmp_${dir}/ ; sed -i "s/BASH_SUB/${Name}$dir/g" pbsscript ; cp   ${sourcedir}/singlepoint.fdf  ${sourcedir}/tmp_${dir}/ ; cat  ${sourcedir}/${dir} >>  ${sourcedir}/tmp_${dir}/singlepoint.fdf ; echo '%endblock' >>  ${sourcedir}/tmp_${dir}/singlepoint.fdf ; cd  ${sourcedir}  ;  done

sourcedir=`pwd` ; for dir in $(ls Plus_* Minus_* ) ; do cd ${sourcedir}/tmp_${dir} ; qsub pbsscript ; done





#once done

list=`ls | grep tmp_Plus` ; for dir in $list ; do tmpvar=`echo ${dir} |  cut -d "_" -f3` ; echo ${tmpvar} ; cp ./${dir}/scf.out  ./plus.out${tmpvar} ; done 

cp plus.out0 minus.out0

list=`ls | grep tmp_Minus` ; for dir in $list ; do tmpvar=`echo ${dir} |  cut -d "_" -f3` ; echo ${tmpvar} ; cp ./${dir}/scf.out  ./minus.out${tmpvar} ; done 



#      DYNMAT


# if Config = 3* Natoms 
/home/pdarancet/SRC/ScriptsSiesta/Raman/extract.sh 28 84

cp forces_plus cut_plus
cp forces_minus cut_minus

# 

cp ../../Si2_EField_500/Phonons/vwhole_junct_NMD.f  ./
#cp ../../Si4_Charged_1000/Phonons/vwhole_junct_NMD.f ./


ifort vwhole_junct_NMD.f -o vwhole_junct_NMD.exe -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread ;   ./vwhole_junct_NMD.exe  4  28 28 84

tail -n 78 out.dat | awk '{print $2, "1.0"}' > phononsDos.dat




~/SRC/ScriptsSiesta/Raman/linei.exe phononsDos.dat 78 phononstmp 3200
~/SRC/ScriptsSiesta/Raman/smear.exe phononstmp 3199 phononsDos.Smeared.dat 10.0













#once done

list=`ls | grep tmp_Plus` ; for dir in $list ; do tmpvar=`echo ${dir} |  cut -d "_" -f3` ; echo ${tmpvar} ; cp ./${dir}/scf.out  ./plus.out${tmpvar} ; done 

cp plus.out0 minus.out0

list=`ls | grep tmp_Minus` ; for dir in $list ; do tmpvar=`echo ${dir} |  cut -d "_" -f3` ; echo ${tmpvar} ; cp ./${dir}/scf.out  ./minus.out${tmpvar} ; done 



#      DYNMAT


# if Config = 3* Natoms 
/home/pdarancet/SRC/ScriptsSiesta/Raman/extract.sh 28 84

cp forces_plus cut_plus
cp forces_minus cut_minus

# if not 
#/home/pdarancet/SRC/ScriptsSiesta/Raman/extract.sh ?? 84 > cut_plus
#cp forces_minus cut_minus



# 

cp ../../Si2_EField_500/Phonons/vwhole_junct_NMD.f  ./
#cp ../../Si4_Charged_1000/Phonons/vwhole_junct_NMD.f ./


ifort vwhole_junct_NMD.f -o vwhole_junct_NMD.exe -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread ;   ./vwhole_junct_NMD.exe  4  28 28 84

tail -n 78 out.dat | awk '{print $2, "1.0"}' > phononsDos.dat




~/SRC/ScriptsSiesta/Raman/linei.exe phononsDos.dat 78 phononstmp 3200
~/SRC/ScriptsSiesta/Raman/smear.exe phononstmp 3199 phononsDos.Smeared.dat 10.0





















#once done

list=`ls | grep tmp_Plus` ; for dir in $list ; do tmpvar=`echo ${dir} |  cut -d "_" -f3` ; echo ${tmpvar} ; cp ./${dir}/scf.out  ./plus.out${tmpvar} ; done 

cp plus.out0 minus.out0

list=`ls | grep tmp_Minus` ; for dir in $list ; do tmpvar=`echo ${dir} |  cut -d "_" -f3` ; echo ${tmpvar} ; cp ./${dir}/scf.out  ./minus.out${tmpvar} ; done 



#      DYNMAT


# if Config = 3* Natoms 
/home/pdarancet/SRC/ScriptsSiesta/Raman/extract.sh 46 138

cp forces_plus cut_plus
cp forces_minus cut_minus


# 



# 

cp ../../Si4_Charged_1000/Phonons/vwhole_junct_NMD.f ./


ifort vwhole_junct_NMD.f -o vwhole_junct_NMD.exe -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread ;   ./vwhole_junct_NMD.exe  4  46 46 138




tail -n 132 out.dat | awk '{print $2, "1.0"}' > phononsDos.dat


~/SRC/ScriptsSiesta/Raman/linei.exe phononsDos.dat 132 phononstmp 3200
~/SRC/ScriptsSiesta/Raman/smear.exe phononstmp 3199 phononsDos.Smeared.dat 10.0

