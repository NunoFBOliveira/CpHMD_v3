#!/bin/bash -e
#
# Read arguments and make some assignments:
[ $# -gt 1 ] && echo "Usage: $0 run_CpHMD.sh" >&2 && exit 1
if [ ! -f $1 ]; then echo "File $1 is missing!!!... Program will crash"; exit 1;
else
    source $1
fi
rundir=`pwd -P`
runfile=`basename $0`
#
# Write executable script:
cat <<EOF > $SysName.slurm
#!/bin/bash -e
InitDate=\`date\`
# Finds current block
export rdir=$rundir
i=1
export j=001

while [ -f ${SysName}_\${j}.occ ] ; do
        i=\$((i+1))
        export j=\`printf "%03d\n" \${i}\`
done

k=\$((i-1))
l=\`printf "%03d\n" \${k}\`

# Info about the Job:
echo "Job executed in Machine \$HOSTNAME" > $rundir/${SysName}_\${j}.blockinfo
echo "Job executed with $nCPU processors" >> $rundir/${SysName}_\${j}.blockinfo
echo "Job executed in DIR: /tmp/${USER}_CpHMD\$$ " >> $rundir/${SysName}_\${j}.blockinfo

echo "" >> $rundir/${SysName}_\${j}.blockinfo
echo -e "Job started on: " >> $rundir/${SysName}_\${j}.blockinfo
date >> $rundir/${SysName}_\${j}.blockinfo

# Copy important files for local directory; first gro is called 
mkdir -p /tmp/${USER}_CpHMD\$$

# Enter local directory
cd /tmp/${USER}_CpHMD\$$

if [ \${i} -eq 1 ]; then 
   cp -f $rundir/$1 ${SysName}_\${j}.pHmdp
elif [[ \${plumedtype} == "grid" ]] ;then
   cp -f $rundir/${SysName}_\${l}.gro ./
   cp -f $rundir/${hills}_\${j} ./
   cp -f $rundir/${colvar_name} ./
   cp -f $rundir/${grid_name} ./
   sed "s/GROin=.*/GROin=${SysName}_\${l}.gro/g" $rundir/$1 > ${SysName}_\${j}.pHmdp
elif [[ \${plumedtype} == "hill" ]] || [[ \${plumedtype} == "static" ]] ;then
   cp -f $rundir/${SysName}_\${l}.gro ./
   cp -f $rundir/${hills} ./
   cp -f $rundir/${colvar_name} ./
   sed "s/GROin=.*/GROin=${SysName}_\${l}.gro/g" $rundir/$1 > ${SysName}_\${j}.pHmdp
else
   cp -f $rundir/${SysName}_\${l}.gro ./
   sed "s/GROin=.*/GROin=${SysName}_\${l}.gro/g" $rundir/$1 > ${SysName}_\${j}.pHmdp
fi

Seg_tot_time=\`awk -v n=${Seg_size} '\$2 ~ /^EffectiveSteps=/ {m=substr (\$2,16); print (n*1000)/(m*0.002)}' ${SysName}_\${j}.pHmdp\`

sed -i "s/InitCycle=.*/InitCycle=\$((\$Seg_tot_time*\$k+1))/g" ${SysName}_\${j}.pHmdp
sed -i "s/EndCycle=.*/EndCycle=\$((\$Seg_tot_time*\$i))/g"  ${SysName}_\${j}.pHmdp

# line needed to mount /programs using autofs
ls /programs/ >/dev/null 2>&1; sleep 2
#
CpHMDDIR=\`egrep "CpHDIR="  ${SysName}_\${j}.pHmdp 2>/dev/null| awk '{print substr(\$2,9,length(\$2)-9)}' 2>/dev/null\`
cp -f \$CpHMDDIR/scripts/CpHMD.sh .

#####################################
## slurm addition to deal with gpu ##
#####################################

if [ $GPU -eq 1 ] ;then	
   OMP_NUM_THREADS=$nCPU

   ### Check if partition is not GPUD ###
   if [[ \$Partition == *"GPUD"* ]] ; then
      echo "ERROR: CpHMD is not suited for running on double GPU setups, please do not use partition GPUD" >> $rundir/${SysName}_\${j}.blockinfo
      exit 1
   fi 

   if [[ \`nvidia-smi -L | wc -l\` > 1 ]] ; then
        if [[ \`ps -aux | awk '{for(i=1;i<=NF;i++){ if(\$i=="-gpu_id"){print 1-\$(i+1)}}}' | wc -l\` > 1 ]]
        then
                echo "Error: no GPU is available for GPU job" >> $rundir/${SysName}_\${j}.blockinfo
                exit 1
        elif [[ \`ps -aux | awk '{for(i=1;i<=NF;i++){ if(\$i=="-gpu_id"){print 1-\$(i+1)}}}' | wc -l\` < 1 ]]
        then
	        export GPUid=0
	else
                export GPUid=\`ps -aux | awk '{for(i=1;i<=NF;i++){ if(\$i=="-gpu_id"){print 1-\$(i+1)}}}'\`
        fi

    else
        export GPUid=0
    fi

    ## Edit the mdrun line in the .pHmdp ##

    sed -i 's/export mdrungpu/export mdrun/' ${SysName}_\${j}.pHmdp
    
else
    sed -i 's/export mdruncpu/export mdrun/' ${SysName}_\${j}.pHmdp
fi

# Run Constant-pH MD segment:
nice -n 19 /tmp/${USER}_CpHMD\$$/CpHMD.sh  ${SysName}_\${j}.pHmdp > ${SysName}_\${j}.err 2>&1 

# Copy files and return to shared directory
cp -df ${SysName}_\${j}* ${rundir}
if [ ! -f ${rundir}/${SysName}-all.sites ] ; then
cp ./${SysName}-all.sites ${rundir}/${SysName}.sites
fi

if [ $ReduceTitration -eq 1 ] ; then
if [ -f ./${SysName}-reducedtitration.sites ] ; then
cp ./${SysName}-reducedtitration.sites ${rundir}/${SysName}_\${j}-RT-debug.dat
fi
if [ -f ./${SysName}.pocc_RT ] ; then
cp ./${SysName}.pocc_RT ${rundir}/${SysName}_\${j}-debug.pocc_RT
fi
fi

if [ -f ./${colvar_name} ] && [ -f ./${hills} ] && [[ \$plumed != "grid" ]]  ; then
cp ./${colvar_name} $rundir
cp ./${hills} $rundir
fi
#
if [ -f ./${colvar_name} ] && [[ \$plumed == "grid" ]]  ; then
cp ./${colvar_name} $rundir
mv ${hills}_curr_seg ${hills}_\${j}
cp ./${hills}_\${j} $rundir
cp ./${grid_name} $rundir 
fi


if (for f in ${SysName}_\${j}*; do diff \$f ${rundir}/\$f; done); then
cd ${rundir}
rm -rf /tmp/${USER}_CpHMD\$$
gzip -9 ${SysName}_\${j}.{err,log,tpr}
else
echo "Error in file copy... please check local files" >> ${rundir}/$rundir/${SysName}_\${j}.blockinfo
exit 1
fi
#
echo "" >> $rundir/${SysName}_\${j}.blockinfo
echo -e "Job finished on: " >> $rundir/${SysName}_\${j}.blockinfo
date >> $rundir/${SysName}_\${j}.blockinfo
#
# Launch next job before exiting:
if [ \${i} -lt $Segments ] # usually 1 segment == 1 ns
then
    cd $rundir; ./$runfile $1
fi
EOF

chmod +x $SysName.slurm

sbatch -p $Partition -N 1 -n $nCPU --mem=1 -o $SysName.sout -e $SysName.serr $SysName.slurm

echo ""
echo "Job submitted to Partition(s): $Partition with $nCPU Processors"
#
## End of Script
#
exit
