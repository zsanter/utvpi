#!/bin/sh

#This is an example script for executing generic jobs with 
# the use of the command 'qsub <name of this script>'


#These commands set up the Grid Environment for your job.  Words surrounding by a backet ('<','>') should be changed
#Any of the PBS directives can be commented out by placing another pound sign in front
#example
##PBS -N name
#The above line will be skipped by qsub because of the two consecutive # signs 

# Specify job name
#PBS -N profile_${program}_${modify}_${f}

# Specify the resources need for the job
# Walltime is specified as hh:mm:ss (hours:minutes:seconds)
#PBS -l nodes=1:ppn=1,pvmem=54gb,walltime=168:00:00

# Specify when Moab should send e-mails 'ae' below user will
# receive e-mail for any errors with the job and/or upon completion
# If you don't want e-mails just comment out these next two PBS lines
#PBS -m ae

# Specify the e-mail address to receive above mentioned e-mails
#PBS -M zsanter@mix.wvu.edu

# Specify the queue to execute task in. Current options can be found by executing the command qstat -q at the terminal
#PBS -q comm_mmem_week

# Enter your command below with arguments just as if you were going to execute on the command line
# It is generally good practice to issue a 'cd' command into the directory that contains the files
# you want to use or use full path names

cd subramani/utvpi
./auto.sh "${program}" "${modify}" "${f}"