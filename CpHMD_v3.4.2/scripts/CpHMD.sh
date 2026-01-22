#!/bin/bash -e 
#
prog=`basename $0` 
usage="Usage: $prog <MyProtein_001>.pHmdp\n
       <MyProtein_001>.pHmdp  : Constant-pH MD parameter file"
#
# This allows to use extended pattern matching features
# E.g.: $PDBin == *(" ")
shopt -s extglob #SC 14-12-2011
#
message ()
{
    case "$1" in
        E ) shift; echo -e "$prog: Error: $*" >&2 ; echo -e "$prog: Error: $*" >> $SLURM_SUBMIT_DIR/${blockname}.blockinfo ; exit 1;;
        U ) shift; echo -e "$prog: Warning: $*\n$usage" >&2; exit 1;;
        W ) shift; echo -e "$prog: Warning: $*" >&2;;
        * ) message E "Wrong use of 'message' function.";;
    esac
}
#
# Parse arguments
if [ $# != 1 ]; then message U "Wrong number of arguments."; fi
#
# Read Simulation Parameters
source $1
#
# Name of your simulation segment (the output files will
# be generated with this name)
export blockname=${1%.*}
export runname=${blockname%_*}
#
# Check some parameters

if [[ ! -d "$ffDIR"/${ffID}.ff ]] ; then
    message E "ffID = $ffID is not valid since the $ffID.ff is not present in your $ffDIR. Check $1 for either your ffID naming or the ffDIR given."
fi

#
# Call all functions stored in the "functions" file:
functions=$CpHDIR/scripts/functions.sh
if [[ -f $functions ]]; then
    source $functions            
else
    message E "File $functions is missing. Check CpHDIR parameter in $1."
fi

######################################################### 
# add the treatment of CpHMD.settings to mdp and fixgro #
#########################################################

awk '/#mdp# / {print}' $1 | sed 's/#mdp# //g' > ${SysName}.mdp
awk '/#fixgro# / {print}' $1 | sed 's/#fixgro# //g' > ${SysName}.fixgro

##################################################### 
# add the treatment of plumed input when plumed =1  #
#####################################################
if  [[ -n $plumed && "${plumed}" -eq "1" ]] ; then
    awk '/#plumed# / {print}' $1 | sed 's/#plumed# //g' > ${SysName}_plumed.dat
    case $plumedtype in
	grid)
	    sed -i "s/\&colvar_stride/$colvar_stride/g" ${SysName}_plumed.dat
	    sed -i "s/\&grid_name/$grid_name/g" ${SysName}_plumed.dat
	    sed -i "s/\&hills/$hills/g" ${SysName}_plumed.dat
	    sed -i "s/\&colvar_name/$colvar_name/g" ${SysName}_plumed.dat
	    ;;
	hill|static)
	    sed -i "s/\&colvar_stride/$colvar_stride/g" ${SysName}_plumed.dat
	    sed -i "s/\&hills/$hills/g" ${SysName}_plumed.dat
	    sed -i "s/\&colvar_name/$colvar_name/g" ${SysName}_plumed.dat
	    ;;
    esac
fi

#
correct_variables #SC 17/3/2010
#
# Let's keep track of the simulation place and time:
echo -e "Simulation run by $USER  @  $HOSTNAME\nInitial time: `date`" \
    > ${blockname}.info
#
# Do some housekeeping
clean_up
#
# Build ff links needed for simulations
build_forcefield
#
# Check if all files are present:
check_files
#
# Make auxiliary files
make_auxiliary_files
#
#Make sites file on the fly
make_sites $1
#
# Get .st files
"$DelphiDir"/getst ${runname}.sites "$StDIR"
#
# Make the delphi database to use
make_delphi_DB
#
#### Starts the constant-pH MD cycle ####
#
TimeStep=`awk '/dt *=/{print $3}' ${runname}.mdp`
#
if [[ $TimeStep == "" ]]; then TimeStep=0.001; fi #SC 25-11-2011
WriteStep=`awk '/nstxout-compressed *=/{print $3}' ${runname}.mdp`
WriteTime=`echo "$WriteStep*$TimeStep" | bc -l`
#
for (( Cycle=$InitCycle ; Cycle <=$EndCycle ; Cycle++ )); do
    #
    # This section is just an index to keep track of the simulation
    sim_time=`echo "$EffectiveSteps*$TimeStep*($Cycle-1)+$WriteTime" | bc -l`
    echo -e "\nCycle = $Cycle; time = $sim_time ps; Date: `date "+%D %T"`" \
         >> ${blockname}.info
    #
    ################### PB/MC PART #####################
    #
    # Make sure the sites file is not empty...
    sitenumball=$(($(wc -l < ${runname}-all.sites)))
    if [ $sitenumball -ne 0 ]; then
	if [ $ReduceTitration == 1 ]; then
	    if [ $((Cycle % RTInterval)) -eq 1 ]; then
		### Get rid of previous reduced sites ###
		rm -f ${runname}.sites
		### make sure the -all sites is now the new one ###
		cp -f ${runname}-all.sites ${runname}.sites
		#####################################
		# ... and call the PB/MC function,...
		echo -n "PB/MC (All) -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
		     >> ${blockname}.info
		run_PBMC red
		echo "`date "+%D %T"`" >> ${blockname}.info
		write_fractions_all_sites
		#
		# ...write fractions to files and build a new topology.
	    else
		sitenumb=$(($(wc -l < ${runname}.sites)))
		if [ $sitenumb -ne 0 ]; then
		    # ... and call the PB/MC function,...
		    echo -n "PB/MC -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
			 >> ${blockname}.info
		    run_PBMC
		    echo "`date "+%D %T"`" >> ${blockname}.info
		    write_fractions
		    
		fi
	    fi
	else
	    sitenumb=$(($(wc -l < ${runname}-all.sites)))
	    if [ $sitenumb -ne 0 ]; then
		# ... and call the PB/MC function,...
		echo -n "PB/MC -        Cycle = $Cycle; Date: `date "+%D %T"` - " \
		     >> ${blockname}.info
		run_PBMC
		echo "`date "+%D %T"`" >> ${blockname}.info
		write_fractions_all_sites
	    fi
	fi
	
	update_topology

	else
        # ... skip the PB/MC and write fractions to files 
        message W "File ${runname}.sites is empty. PB/MC step is not performed in cycle $Cycle."
        echo "" >> "TMP_CpHMD.occ"
        echo "" >> "TMP_CpHMD.mocc"
	#migrate the charge topology back to its original name
	mv TMP_processed.top TMP_CpHMD.top
    fi
    
    #### MD PART ####
    # Call dynamics with solvent relaxation 
    if [ $RelaxSteps != 0 ]; then
        echo -n "MD relax     - Cycle = $Cycle; Date: `date "+%D %T"` - " \
            >> ${blockname}.info
	run_relaxation #SC 28-11-2011
        echo "`date "+%D %T"`" >> ${blockname}.info
    else
        mv TMP_effective.gro TMP_relax.gro
    fi
    #
    # Call effective (full) dynamics
    echo -n "MD effective - Cycle = $Cycle; Date: `date "+%D %T"` - " \
        >> ${blockname}.info
    run_dynamics effective relax
    echo "`date "+%D %T"`" >> ${blockname}.info
    #
         

    # Call Append data function   
    data_append
done
#### Ends the constant-pH MD cycle ####
#
# Store Segment Outputs with unambigous Name
for e in gro tpr edr log xtc; do  mv -f TMP_CpHMD.$e ${blockname}.$e; done
## move the energy calculation
if [ -f Eb_calculation.dat ]; then
    mv Eb_calculation.dat ${blockname}.ene
fi
#
if [ -f TMP_CpHMD.occ ]; then
    for e in occ mocc ; do mv -f TMP_CpHMD.$e ${blockname}.$e; done
fi

if [ -f TMP_CpHMD_pullx.xvg ]; then
    for e in x f ; do mv -f TMP_CpHMD_pull$e.xvg ${blockname}_pull$e.xvg; done
fi
# Correct final timestamps of hills and colvar:
if [ -f $colvar_name ]; then
    sed -i '/^ 20.0/d' ${colvar_name}  
fi
#
# Let's keep track of the simulation end time:
echo -e "\nEnd time:     `date`" >> ${blockname}.info
#
# Clean up function

clean_up

#
exit 0
