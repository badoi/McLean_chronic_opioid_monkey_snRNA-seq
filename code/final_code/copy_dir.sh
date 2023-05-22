PROJDIR=/projects/pfenninggroup/singleCell/McLean_chronic_opioid_monkey_snRNA-seq
PROJDIR2=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq

rsync -Pav ${PROJDIR2}/code/final_code/HeKleyman2021_macaque_striatum_data_processing \
${PROJDIR}/code/final_code

rsync -Pav ${PROJDIR2}/data/tidy_data/HeKleyman2021_macaque_striatum_data_processing \
${PROJDIR}/data/tidy_data


rsync -Pav ${PROJDIR2}/code/raw_code/preprocess_snRNA-seq_reads \
${PROJDIR}/code/raw_code

