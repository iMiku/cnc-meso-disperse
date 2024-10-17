#!/bin/bash

# List of strings
vol_fracs=("0.0025" "0.0050" "0.0075" "0.0100" "0.0125" "0.0150")
aspect_ratios=("010" "030" "050" "070" "090")
charges=("0.00" "0.25" "0.50" "0.75" "1.00" "1.25" "1.50" "1.75")
fstds=("30.0")
gn=6
lmp_cmd_pre="mpirun -n 8 path/lmp_mpi"
lmp_cmd_suf="-screen none -in cnc_disperse_wat.in"
# Loop over the list
for vol_frac in "${vol_fracs[@]}"; do
	for aspect_ratio in "${aspect_ratios[@]}"; do
		for charge in "${charges[@]}"; do
			for fstd in "${fstds[@]}"; do
				start_time=$(date +%s)

				cmd_para1="-var gn ${gn} -var fstd ${fstd} -var charge ${charge}"
				cmd_para2="-var vol_frac ${vol_frac} -var aspect_ratio ${aspect_ratio}"
				cmd="${lmp_cmd_pre} ${cmd_para1} ${cmd_para2} ${lmp_cmd_suf}"
				echo ${cmd}
				${cmd}

				end_time=$(date +%s)
				duration=$((end_time - start_time))
				echo "Time taken for ${duration} seconds"
			done
		done
	done
done
