#!/bin/bash 

preis="/home/noliveira/08_USCpHMD-preis/pHRE/pentapeptides/"

echo -n -e "Which residues to check? (all | asp glu his lys cys tyr cter nter)\n"
read residues

if [ $residues == "all" ]
then
    residues="asp glu his lys cys tyr cter nter"
fi

for res in $residues  #asp glu his lys cys tyr cter nter
do
    #pH_list=`ls -a $preis/${res}_4reps_20ps_r1 | grep -E pH[0-9] `
    pH_list=`ls -a ./${res} | grep -E pH[0-9] `
    echo "-----------------------------------------------------"
    for pH in pH05.00 #$pH_list
    do		
	for i in 2 3 4 5
	do
	    echo " "
	    echo -n -e "Residue - ${res} ${pH} rep${i}\n\n"
	    
	    if [ -d "./${res}/${pH}/rep${i}" ]
	    then
		
		ls -a ./${res}/${pH}/rep${i} | grep d_ > ./tmp.txt
		
		div=`cat tmp.txt | wc -l | awk '{print int($1/4)}'`

		#echo $div
		for ((l=1;l<=$div+1;l++))
		do
		    
		    first=`awk -v l=$l  'NR==((l*4)-3) {print}' tmp.txt`
		    second=`awk -v l=$l 'NR==((l*4)-2) {print}' tmp.txt`
		    third=`awk -v l=$l  'NR==((l*4)-1) {print}' tmp.txt`
		    fourth=`awk -v l=$l 'NR==((l*4)) {print}' tmp.txt`
		    #fifth=`awk -v l=$l 'NR==((l*5)) {print}' tmp.txt`

		    #echo $first $second $third
		    
		    for d in $first $second $third $fourth
		    do
			
			sim=`s-id noliveira | grep ${res}/${pH}/rep${i}/${d}`

			if [ -z "$sim" ] 
			then		    		    
        		    nano=`ls -1 ./${res}/${pH}/rep$i/${d}/*.xtc 2>/dev/null | wc -l | awk '{printf("%.0f", $1)}'`
			    
			    if [ $nano -eq "100" ]
			    then
				status=`echo -e "Finish" | awk '{printf("%-8s", "\033[4;36m"$1"\033[0m")}'`
				ns=`ls -1 ./${res}/${pH}/rep$i/${d}/*.xtc 2>/dev/null | wc -l | awk '{printf("%.0f", $1)}'`
				
			    else
				
				status=`echo -e "Stopped" | awk '{printf("%-9s", "\033[31m"$1"\033[0m")}'`
				ns=`ls -1 ./${res}/${pH}/rep$i/${d}/*.xtc 2>/dev/null | wc -l | awk '{printf("%.0f", $1)}'`
			    fi
			    
			    
			    
			else
			    run_Check=`echo $sim | awk '{print $3}'`
			    f=0
			    e=0
			    c=0			
			    if [ $run_Check == "R" ]
			    then			    
				status=`echo -e "Running " | awk '{printf("%-9s", "\033[32m"$1"\033[0m")}'`
				ns=`ls -1 ./${res}/${pH}/rep$i/${d}/*.xtc 2>/dev/null | wc -l | awk '{printf("%.0f", $1)}'`
				
			    elif [ $run_Check == "PD" ] 
			    then
				
				status=`echo -e "Waiting " | awk '{printf("%-9s", "\033[93m"$1"\033[0m")}'`
				ns=`ls -1 ./${res}/${pH}/rep$i/${d}/*.xtc 2>/dev/null | wc -l | awk '{printf("%.0f", $1)}'`
			    else
				status=`echo -e "??? " | awk '{printf("%-9s", $1)}'`
				ns=`ls -1 ./${res}/${pH}/rep$i/${d}/*.xtc 2>/dev/null | wc -l | awk '{printf("%.0f", $1)}'`
				
			    fi
			    
			    
			    #echo $status 
			    
			    
			fi
			echo -n -e "US-$d - $status -=-> $ns \t"
		    done
		    echo ""
		done
	    else
		#echo -n -e "Residue - ${res} ${pH} not started\n\n"
		#echo "Not started"
		echo ""
	    fi
	    
	done
    done
done

