Assess the impact of population structure on the estimation of the rate of adaptation. Simulations run by Prof. Julien Yann Dutheil (https://www.evolbio.mpg.de/person/40579/15929)
	

```bash
NTHREADS=20
```

Panmicic population
===================

Write a SLURM bash script:
```bash
mkdir -p log/err
mkdir -p log/out
SCRIPT=slurm_slim_panmictic.sh
rm $SCRIPT
echo "#! /bin/bash"                                   >> $SCRIPT
echo ""                                               >> $SCRIPT
echo "#SBATCH --job-name=slim_panmictic"              >> $SCRIPT
echo "#SBATCH --ntasks=1"                             >> $SCRIPT
echo "#SBATCH --nodes=1"                              >> $SCRIPT
echo "#SBATCH --array=1-50"                           >> $SCRIPT
echo "#SBATCH --time=2-00:00:00"                      >> $SCRIPT
echo "#SBATCH --mem=64G"                              >> $SCRIPT
echo "#SBATCH --error=log/err/slim_panmictic.%J.err"  >> $SCRIPT
echo "#SBATCH --output=log/out/slim_panmictic.%J.out" >> $SCRIPT
echo "#SBATCH --mail-type=ALL"                        >> $SCRIPT
echo "#SBATCH --mail-user=dutheil@evolbio.mpg.de"     >> $SCRIPT
echo "#SBATCH --partition=medium"                     >> $SCRIPT
echo ""                                               >> $SCRIPT
allPs=(0 0.00001 0.00002 0.0001 0.001)
JOBCOUNT=0
for j in ${!allPs[@]}; do 
  for i in {1..10}; do
    JOBCOUNT=$((JOBCOUNT+1))
    echo "if [ \$SLURM_ARRAY_TASK_ID == $JOBCOUNT ]; then" >> $SCRIPT
    echo "  mkdir -p Output/Panmictic/P$j"      >> $SCRIPT
    echo "  ./slim -seed $((42+i)) -d propPos=${allPs[$j]} SimulatePanmictic.slim > Output/Panmictic/P$j/sim_p${j}_rep${i}.slim.out" >> $SCRIPT
    echo "fi" >> $SCRIPT
  done
done
```

Parse the simulation output:
```bash
allPs=(0 0.00001 0.00002 0.0001 0.001)
for j in ${!allPs[@]}; do 
  for i in {1..10}; do
    ./slim2dfem/slim2dfem Output/Panmictic/P$j/sim_p${j}_rep${i}.slim.out Output/Panmictic/P$j/sim_p${j}_rep${i}_unfolded.dofe unfolded >& Output/Panmictic/P$j/sim_p${j}_rep${i}_unfolded.log
  done
done
```

Now run grapes:
```bash
rm cmd_grapes_panmictic.sh
allPs=(0 0.00001 0.00002 0.0001 0.001)
for j in ${!allPs[@]}; do 
  for i in {1..10}; do
    #CMD="export DYLD_LIBRARY_PATH=$HOME/.local/lib; grapes " #on mac
    CMD="grapes" #on linux
    CMD="$CMD -in Output/Panmictic/P$j/sim_p${j}_rep${i}_unfolded.dofe"
    CMD="$CMD -out Output/Panmictic/P$j/sim_p${j}_rep${i}_unfolded_grapes.csv"
    CMD="$CMD -model GammaExpo -no_div_param -nb_rand_start 5 -nearly_neutral 0."
    CMD="$CMD > Output/Panmictic/P$j/sim_p${j}_rep${i}_unfolded_grapes.log"
    echo $CMD >> cmd_grapes_panmictic.sh
  done
done
parallel -j $NTHREADS --eta < cmd_grapes_panmictic.sh
```

Population exponential growth/decline
=====================================


Write a SLURM bash script:
```bash
mkdir -p log/err
mkdir -p log/out
SCRIPT=slurm_slim_exponential.sh
rm $SCRIPT
echo "#! /bin/bash"                                     >> $SCRIPT
echo ""                                                 >> $SCRIPT
echo "#SBATCH --job-name=slim_exponential"              >> $SCRIPT
echo "#SBATCH --ntasks=1"                               >> $SCRIPT
echo "#SBATCH --nodes=1"                                >> $SCRIPT
echo "#SBATCH --array=1-250"                            >> $SCRIPT
echo "#SBATCH --time=2-00:00:00"                        >> $SCRIPT
echo "#SBATCH --mem=64G"                                >> $SCRIPT
echo "#SBATCH --error=log/err/slim_exponential.%J.err"  >> $SCRIPT
echo "#SBATCH --output=log/out/slim_exponential.%J.out" >> $SCRIPT
echo "#SBATCH --mail-type=ALL"                          >> $SCRIPT
echo "#SBATCH --mail-user=dutheil@evolbio.mpg.de"       >> $SCRIPT
echo "#SBATCH --partition=medium"                       >> $SCRIPT
echo ""                                                 >> $SCRIPT
allPs=(0 0.00001 0.00002 0.0001 0.001)
allNs=(1000 5000 10000 20000 100000)
JOBCOUNT=0
for k in ${!allNs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      JOBCOUNT=$((JOBCOUNT+1))
      echo "if [ \$SLURM_ARRAY_TASK_ID == $JOBCOUNT ]; then" >> $SCRIPT
      echo "  mkdir -p Output/Exponential/N$k/P$j"      >> $SCRIPT
      echo "  ./slim -seed $((42+i)) -d propPos=${allPs[$j]} -d newPopSize=${allNs[$k]} SimulateExponential.slim > Output/Exponential/N$k/P$j/sim_n${k}_p${j}_rep${i}.slim.out" >> $SCRIPT
      echo "fi" >> $SCRIPT
    done
  done
done
```

Parse the simulation output:
```bash
allPs=(0 0.00001 0.00002 0.0001 0.001)
allNs=(1000 5000 10000 20000 100000)
for k in ${!allNs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      ./slim2dfem/slim2dfem Output/Exponential/N$k/P$j/sim_n${k}_p${j}_rep${i}.slim.out \
                            Output/Exponential/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded.dofe \
                unfolded >& Output/Exponential/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded.log
    done
  done
done
```

Now run grapes:
```bash
rm cmd_grapes_exponential.sh
allPs=(0 0.00001 0.00002 0.0001 0.001)
allNs=(1000 5000 10000 20000 100000)
for k in ${!allNs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      #CMD="export DYLD_LIBRARY_PATH=$HOME/.local/lib; grapes " #on mac
      CMD="grapes" #on linux
      CMD="$CMD -in Output/Exponential/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded.dofe"
      CMD="$CMD -out Output/Exponential/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded_grapes.csv"
      CMD="$CMD -model GammaExpo -no_div_param -nb_rand_start 5 -nearly_neutral 0."
      CMD="$CMD > Output/Exponential/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded_grapes.log"
      echo $CMD >> cmd_grapes_exponential.sh
    done
  done
done
parallel -j $NTHREADS --eta < cmd_grapes_exponential.sh
```


Simple population structure: two demes, constant
================================================

Write a SLURM bash script:
```bash
mkdir -p log/err
mkdir -p log/out
SCRIPT=slurm_slim_2demes.sh
rm $SCRIPT
echo "#! /bin/bash"                                >> $SCRIPT
echo ""                                            >> $SCRIPT
echo "#SBATCH --job-name=slim_2demes"              >> $SCRIPT
echo "#SBATCH --ntasks=1"                          >> $SCRIPT
echo "#SBATCH --nodes=1"                           >> $SCRIPT
echo "#SBATCH --array=1-250"                       >> $SCRIPT
echo "#SBATCH --time=2-00:00:00"                   >> $SCRIPT
echo "#SBATCH --mem=64G"                           >> $SCRIPT
echo "#SBATCH --error=log/err/slim_2demes.%J.err"  >> $SCRIPT
echo "#SBATCH --output=log/out/slim_2demes.%J.out" >> $SCRIPT
echo "#SBATCH --mail-type=ALL"                     >> $SCRIPT
echo "#SBATCH --mail-user=dutheil@evolbio.mpg.de"  >> $SCRIPT
echo "#SBATCH --partition=medium"                  >> $SCRIPT
echo ""                                            >> $SCRIPT
allPs=(0 0.00001 0.00002 0.0001 0.001)
allAs=(0.0000001 0.000001 0.00001 0.0001 0.001)
JOBCOUNT=0
for k in ${!allAs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      JOBCOUNT=$((JOBCOUNT+1))
      echo "if [ \$SLURM_ARRAY_TASK_ID == $JOBCOUNT ]; then" >> $SCRIPT
      echo "  mkdir -p Output/TwoDemes/A$k/P$j"      >> $SCRIPT
      echo "  ./slim -seed $((42+i)) -d propPos=${allPs[$j]} -d propAdmixture=${allAs[$k]} SimulateTwoDemes.slim > Output/TwoDemes/A$k/P$j/sim_a${k}_p${j}_rep${i}.slim.out" >> $SCRIPT
      echo "fi" >> $SCRIPT
    done
  done
done
```

Parse the simulation output:
```bash
allPs=(0 0.00001 0.00002 0.0001 0.001)
allAs=(0.0000001 0.000001 0.00001 0.0001 0.001)
for k in ${!allAs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      ./slim2dfem/slim2dfem Output/TwoDemes/A$k/P$j/sim_a${k}_p${j}_rep${i}.slim.out \
                            Output/TwoDemes/A$k/P$j/sim_a${k}_p${j}_rep${i}_unfolded.dofe \
                unfolded >& Output/TwoDemes/A$k/P$j/sim_a${k}_p${j}_rep${i}_unfolded.log
    done
  done
done
```

Now run grapes:
```bash
rm cmd_grapes_2demes.sh
allPs=(0 0.00001 0.00002 0.0001 0.001)
allAs=(0.0000001 0.000001 0.00001 0.0001 0.001)
for k in ${!allAs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      #CMD="export DYLD_LIBRARY_PATH=$HOME/.local/lib; grapes " #on mac
      CMD="grapes" #on linux
      CMD="$CMD -in Output/TwoDemes/A$k/P$j/sim_a${k}_p${j}_rep${i}_unfolded.dofe"
      CMD="$CMD -out Output/TwoDemes/A$k/P$j/sim_a${k}_p${j}_rep${i}_unfolded_grapes.csv"
      CMD="$CMD -model GammaExpo -no_div_param -nb_rand_start 5 -nearly_neutral 0."
      CMD="$CMD > Output/TwoDemes/A$k/P$j/sim_a${k}_p${j}_rep${i}_unfolded_grapes.log"
      echo $CMD >> cmd_grapes_2demes.sh
    done
  done
done
parallel -j $NTHREADS --eta < cmd_grapes_2demes.sh
```

Population structure: 2 demes merge
===================================

Write a SLURM bash script:
```bash
mkdir -p log/err
mkdir -p log/out
SCRIPT=slurm_slim_2demes_merge.sh
rm $SCRIPT
echo "#! /bin/bash"                                      >> $SCRIPT
echo ""                                                  >> $SCRIPT
echo "#SBATCH --job-name=slim_2demes_merge"              >> $SCRIPT
echo "#SBATCH --ntasks=1"                                >> $SCRIPT
echo "#SBATCH --nodes=1"                                 >> $SCRIPT
echo "#SBATCH --array=1-250"                             >> $SCRIPT
echo "#SBATCH --time=2-00:00:00"                         >> $SCRIPT
echo "#SBATCH --mem=64G"                                 >> $SCRIPT
echo "#SBATCH --error=log/err/slim_2demes_merge.%J.err"  >> $SCRIPT
echo "#SBATCH --output=log/out/slim_2demes_merge.%J.out" >> $SCRIPT
echo "#SBATCH --mail-type=ALL"                           >> $SCRIPT
echo "#SBATCH --mail-user=dutheil@evolbio.mpg.de"        >> $SCRIPT
echo "#SBATCH --partition=medium"                        >> $SCRIPT
echo ""                                                  >> $SCRIPT
allPs=(0 0.00001 0.00002 0.0001 0.001)
allMs=(20000 40000 60000 80000 100000)
JOBCOUNT=0
for k in ${!allMs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      JOBCOUNT=$((JOBCOUNT+1))
      echo "if [ \$SLURM_ARRAY_TASK_ID == $JOBCOUNT ]; then" >> $SCRIPT
      echo "  mkdir -p Output/TwoDemesMerge/M$k/P$j"         >> $SCRIPT
      echo "  ./slim -seed $((42+i)) -d propPos=${allPs[$j]} -d propAdmixture=0.000001 SimulateTwoDemesMerge${allMs[$k]}.slim > Output/TwoDemesMerge/M$k/P$j/sim_m${k}_p${j}_rep${i}.slim.out" >> $SCRIPT
      echo "fi" >> $SCRIPT
    done
  done
done
```

Parse the simulation output:
```bash
allPs=(0 0.00001 0.00002 0.0001 0.001)
allMs=(20000 40000 60000 80000 100000)
for k in ${!allAs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      ./slim2dfem/slim2dfem Output/TwoDemesMerge/M$k/P$j/sim_m${k}_p${j}_rep${i}.slim.out \
                            Output/TwoDemesMerge/M$k/P$j/sim_m${k}_p${j}_rep${i}_unfolded.dofe \
                unfolded >& Output/TwoDemesMerge/M$k/P$j/sim_m${k}_p${j}_rep${i}_unfolded.log
    done
  done
done
```

Now run grapes:
```bash
rm cmd_grapes_2demes_merge.sh
allPs=(0 0.00001 0.00002 0.0001 0.001)
allMs=(20000 40000 60000 80000 100000)
for k in ${!allAs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      #CMD="export DYLD_LIBRARY_PATH=$HOME/.local/lib; grapes " #on mac
      CMD="grapes" #on linux
      CMD="$CMD -in Output/TwoDemesMerge/M$k/P$j/sim_m${k}_p${j}_rep${i}_unfolded.dofe"
      CMD="$CMD -out Output/TwoDemesMerge/M$k/P$j/sim_m${k}_p${j}_rep${i}_unfolded_grapes.csv"
      CMD="$CMD -model GammaExpo -no_div_param -nb_rand_start 5 -nearly_neutral 0."
      CMD="$CMD > Output/TwoDemesMerge/M$k/P$j/sim_m${k}_p${j}_rep${i}_unfolded_grapes.log"
      echo $CMD >> cmd_grapes_2demes_merge.sh
    done
  done
done
parallel -j $NTHREADS --eta < cmd_grapes_2demes_merge.sh
```

Population exponential growth/decline with two demes
====================================================


Write a SLURM bash script:
```bash
mkdir -p log/err
mkdir -p log/out
SCRIPT=slurm_slim_exponential_2demes.sh
rm $SCRIPT
echo "#! /bin/bash"                                            >> $SCRIPT
echo ""                                                        >> $SCRIPT
echo "#SBATCH --job-name=slim_exponential_2demes"              >> $SCRIPT
echo "#SBATCH --ntasks=1"                                      >> $SCRIPT
echo "#SBATCH --nodes=1"                                       >> $SCRIPT
echo "#SBATCH --array=1-250"                                   >> $SCRIPT
echo "#SBATCH --time=2-00:00:00"                               >> $SCRIPT
echo "#SBATCH --mem=64G"                                       >> $SCRIPT
echo "#SBATCH --error=log/err/slim_exponential_2demes.%J.err"  >> $SCRIPT
echo "#SBATCH --output=log/out/slim_exponential_2demes.%J.out" >> $SCRIPT
echo "#SBATCH --mail-type=ALL"                                 >> $SCRIPT
echo "#SBATCH --mail-user=dutheil@evolbio.mpg.de"              >> $SCRIPT
echo "#SBATCH --partition=medium"                              >> $SCRIPT
echo ""                                                        >> $SCRIPT
allPs=(0 0.00001 0.00002 0.0001 0.001)
allNs=(1000 5000 10000 20000 100000)
JOBCOUNT=0
for k in ${!allNs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      JOBCOUNT=$((JOBCOUNT+1))
      echo "if [ \$SLURM_ARRAY_TASK_ID == $JOBCOUNT ]; then" >> $SCRIPT
      echo "  mkdir -p Output/ExponentialTwoDemes/N$k/P$j"      >> $SCRIPT
      echo "  ./slim -seed $((42+i)) -d propPos=${allPs[$j]} -d newPopSize=${allNs[$k]} -d propAdmixture=0.000001 SimulateExponentialTwoDemes.slim > Output/ExponentialTwoDemes/N$k/P$j/sim_n${k}_p${j}_rep${i}.slim.out" >> $SCRIPT
      echo "fi" >> $SCRIPT
    done
  done
done
```

Parse the simulation output:
```bash
allPs=(0 0.00001 0.00002 0.0001 0.001)
allNs=(1000 5000 10000 20000 100000)
for k in ${!allNs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      ./slim2dfem/slim2dfem Output/ExponentialTwoDemes/N$k/P$j/sim_n${k}_p${j}_rep${i}.slim.out \
                            Output/ExponentialTwoDemes/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded.dofe \
                unfolded >& Output/ExponentialTwoDemes/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded.log
    done
  done
done
```

Now run grapes:
```bash
rm cmd_grapes_exponential_2demes.sh
allPs=(0 0.00001 0.00002 0.0001 0.001)
allNs=(1000 5000 10000 20000 100000)
for k in ${!allNs[@]}; do 
  for j in ${!allPs[@]}; do 
    for i in {1..10}; do
      #CMD="export DYLD_LIBRARY_PATH=$HOME/.local/lib; grapes " #on mac
      CMD="grapes" #on linux
      CMD="$CMD -in Output/ExponentialTwoDemes/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded.dofe"
      CMD="$CMD -out Output/ExponentialTwoDemes/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded_grapes.csv"
      CMD="$CMD -model GammaExpo -no_div_param -nb_rand_start 5 -nearly_neutral 0."
      CMD="$CMD > Output/ExponentialTwoDemes/N$k/P$j/sim_n${k}_p${j}_rep${i}_unfolded_grapes.log"
      echo $CMD >> cmd_grapes_exponential_2demes.sh
    done
  done
done
parallel -j $NTHREADS --eta < cmd_grapes_exponential_2demes.sh
```
8 runs did not converge after 48h and were killed.



