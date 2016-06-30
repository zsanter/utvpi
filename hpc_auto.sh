#!/bin/sh

# qsub -v program=[program],modify=["variables"|"constraints"],f=["f0"|"f1"|"f2"] hpc_auto.sh

# Specify job name
#PBS -N utvpi_profile

# Specify the resources needed for the job
# 1 node, 1 processor per node, 54 GB of memory, 10 hours max
#PBS -l nodes=1:ppn=1,pvmem=54gb,walltime=10:00:00

# Specify when Moab should send e-mails
# Upon abort/begin/end
#PBS -m abe

# Specify the e-mail address to receive above mentioned e-mails
#PBS -M zsanter@mix.wvu.edu

# Specify the queue to execute task in.
#PBS -q comm_mmem_day

# Specify that the job has exclusive access to the node on which it executes.
#PBS -n

# Whenever a job begins, the node starts in the user's home directory.
cd subramani/utvpi
./auto.sh "${program}" "${modify}" "${f}"