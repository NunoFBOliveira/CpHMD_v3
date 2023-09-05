
#grom=/gromacs/gromacs-2018.6/bin/gmx
grom=/gromacs/gromacs-5.1.5_pH_I/bin/gmx
Dir=`pwd`

zzmin=-1.9
zzmax=0.6

min=10000   
max=250000

zprof0=-4.2
bootstrap=10

/bin/rm -f aux*

for pH in 07.0 #05.0
do
    rm -f tpr-files.dat
    for u in `awk '{print $1}' ../frame_list.dat`
    do
	if [ -f $Dir/pH${pH}/d_${u}/carrier_001.tpr.gz ]; then
	    gunzip $Dir/pH${pH}/d_${u}/carrier_001.tpr.gz
	fi
	ls -1 $Dir/pH${pH}/d_${u}/carrier_001.tpr >> tpr-files.dat

	cat $Dir/pH${pH}/d_${u}/carrier_???_pullf.xvg | \
            awk -v min=$min -v max=$max '$1 > min && $1 <= max' > aux_${u}_f.xvg
        ls -1 aux_${u}_f.xvg >> aux_f-files.dat
    done
    
    cd $Dir
    # Global PMFs
    ${grom} wham -it tpr-files.dat \
	    -if aux_f-files.dat \
	    -o aux_profilef_pH${pH}.xvg \
	    -hist aux_histf_pH${pH}.xvg \
	    -bsres pmf_pH${pH}.xvg \
	    -bsprof pmf_bsprofilesf_pH${pH}.xvg \
	    -unit kCal -xvg none\
	    -tol 1e-6 \
	    -bins 100 \
	    -zprof0 $zprof0 \
	    -b $min \
	    -e $max \
	    -min $zzmin \
	    -max $zzmax \
	    -nBootstrap $bootstrap
    
    rm -f aux_* \#* *\~
done

gpfile='plot.gp'
cat <<EOF >> $gpfile
set term postscript enhanced color solid "Helvetica" 20

set encoding iso_8859_1
set border 31 lt -1 lw 1

set style line  1 lt rgb "#FF0000" lw 3 pt 7 ps 1.0
set style line  2 lt rgb "#0000FF" lw 3 pt 7 ps 1.0
set style line  3 lt rgb "#000000" lw 3 pt 5 ps 1.0

set grid back

set xlabel "Distance do the channel center / nm" font "Helvetica, 25"
set xrange [-4.5:+4.5] 
set xtics -5,1.0
set mxtics 5
set ylabel "Free energy / kcal mol^{-1}" font "Helvetica, 25"
unset yrange
#set yrange [-5:7]
set ytics -120,10
set mytics 2
set key left
set output "plots_profiles.ps"
set xzeroaxis

set object 1 rect from   -2.1,-100 to  2.1,100 fc rgb "#CCCCCC" fillstyle noborder behind

set label "ADP/ATP Transporter" at  -1.8,30 font "Helvetica-Bold,22"
set label "Cytoplasm" at -4.1,+5 font "Helvetica-Bold,22" rotate by 90
set label "Matrix"    at +4.1,-2 font "Helvetica-Bold,22" rotate by 270

plot 'pmf_pH07.0.xvg'  u 1:2   notitle      w l        ls 2 ,\
     'pmf_pH07.0.xvg'  u 1:2:3 every 5 title 'pH 7' w errorbar ls 2 ,\
#     'pmf_pH05.0.xvg'  u 1:2   notitle      w l        ls 1 ,\
#     'pmf_pH05.0.xvg'  u 1:2:3 every 5 title 'pH 5' w errorbar ls 1 ,\


EOF

gnuplot $gpfile
! ps2pdf plots_profiles.ps

rm -f pull_*.xvg *.ps *0_seg-*.xvg *bsprofilesf*  tpr-files.dat $gpfile



