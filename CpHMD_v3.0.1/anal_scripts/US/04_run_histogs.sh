path=`pwd`

mkdir -p Analyses
mkdir -p Analyses/histogs
mkdir -p Analyses/histogs/data

for pH in 07.0
do
    for dist in `awk '{print $1}' ../frame_list-tot.dat`
    do
	echo $dist
	if [ ! -f "$path/pH${pH}/d_${dist}/*pullx.xvg > Analyses/histogs/data/d_${dist}" ]
	then
	    
	    awk '{print NR/10, $2*10}' $path/pH${pH}/d_${dist}/*pullx.xvg > Analyses/histogs/data/d_${dist}

	fi
	
	awk '$1>=20000 {print $2}' Analyses/histogs/data/d_${dist}| histog -r -42.5,42.5 -d 0.1 > $path/Analyses/histogs/histog_${dist}

	tot=`awk '{s+=$2};END{print s}' $path/Analyses/histogs/histog_${dist} `
	echo $tot
	awk -v tot=$tot '{print $1, ($2/tot)/0.1}' $path/Analyses/histogs/histog_${dist} > $path/Analyses/histogs/histog_${dist}_norm
	
	rm -f \#* *\~
    done
done

cd $path/Analyses/histogs/

gpfile=US_histog.gp
psfile=US_histog.ps
    
cat<<EOF > $gpfile
set term postscript enhanced color solid "Helvetica" 22
set output "$psfile"
set tics nomirror

unset key

set pointsize 1


filelist=system("ls -1 ./histog_*_norm")
plot for [data in filelist] data using 1:2 w lines lw 3
EOF
	
gnuplot $gpfile
rm -f aux*
ps2pdf $psfile

#convert -trim -density 150 $file.pdf $file.png

rm $gpfile $psfile

cd $path
