#!/bin/bash


covar_headers=$1
exposures=$2

covar_arr=(${covar_headers})
exp_arr=(${exposures})
new_cov_arr=(${exp_arr[@]})

for cov in ${covar_arr[@]}; do
        is_exp=0
        for exp in ${exp_arr[@]}; do
                if [ ${cov} = ${exp} ]; then is_exp=1; fi
        done
        if [ ${is_exp} = 0 ]; then new_cov_arr+=(${cov}); fi
done

echo ${new_cov_arr[@]}
echo ${#exp_arr[@]}
