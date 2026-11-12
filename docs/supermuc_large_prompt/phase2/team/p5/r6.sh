bash $SCRATCH/large_prompt/team/p5/jobs/t_build.sh 2>&1 | grep -E "spilled|real" | sed -e 's/_ZTSZZN8alfdevqr//' -e 's/EvRN4sycl.*compiled/ compiled/' | cut -c1-140
bash $SCRATCH/large_prompt/team/p5/jobs/t_abl.sh
