#!bin/bash
constrains=( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "24")
n_c=(2 2 3 3 4 4 5 5 5 6 7 7 7 7 8 8 8 8 8 9 9 9 9 10)
#constrains=( "07" "08")
#n_c=(5 5)
#model=('TD')
model=( "d1s1a2" ) #'d4s3a2' 'd7s2a2')
#model=("d1s1" "d4s3" "d7s3")
#model=("d1s1" "d1s3" "d2s1" "d2s2" "d4s3" "d3s2")
for k in ${!model[@]};
do
for j in ${!constrains[@]};
do
log=${constrains[$j]}
log1=${n_c[$j]}
log2=${model[$k]}
#python compare_wmap_planck_ffp10_sim_ILC_dust_pysm.py $log $log1 $log2
#python compare_wmap_planck_ffp10_sim_ILC_dust_pysm_cp.py $log $log1 $log2
python compare_wmap_planck_sim_ILC_dust_pysm.py $log $log1 $log2 'd1' 's1'
#python compare_wmap_planck_sim_ILC_dust_TD.py $log $log1 $log2 'd1' 's1'
#python compare_wmap_planck_sim_ILC_dust_pysm_new_pivot.py $log $log1 $log2 'd1' 's1'
#python compare_wmap_planck_sim_ILC_dust_pysm.py $log $log1 'd4s3a2' 'd4' 's3'
#python compare_wmap_planck_sim_ILC_dust_pysm.py $log $log1 'd7s2a2' 'd7' 's2'


#python compare_wmap_planck_sim_ILC_dust_pysm_v1.py $log $log1 'd1s1' 'd1' 's1'
#python compare_wmap_planck_sim_ILC_dust_pysm_v1.py $log $log1 'd4s3' 'd4' 's3'
#python compare_wmap_planck_sim_ILC_dust_pysm_v1.py $log $log1 'd7s2' 'd7' 's2'
#python compare_wmap_planck_sim_ILC_dust_pysm.py $log $log1 'd7s1' 'd7' 's1'
#python compare_wmap_planck_data_ILC_dust.py $log $log1
done
done
