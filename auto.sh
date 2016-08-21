#!/bin/bash

# ./auto.sh [program] [modify:"variables"|"constraints"] [feasibility:"f0"|"f1"|"f2"]
#
# This script profiles the input program with sets of input systems that either modulate the number of variables, keeping the 
# number of constraints constant; or modulate the number of constraints, keeping the number of variables constant. The feasibility
# input determines the set of input systems, based on the feasibility setting used to generate these systems, which does not 
# correlate directly to the similarly-named "f" output that enters the empirical csv file.
#
# f0 input systems were generated with no guarantee of feasibility; none of the f0 input systems are linearly feasible. f1 input 
# systems were generated with a guarantee of at least linear feasibility; one of these turned out to not be linearly feasible and
# about half of them are also integrally feasible. f2 input systems were generated with a guarantee of integral feasibility; one 
# of these also turned out to not be linearly feasible, though the rest are integrally feasible.
#
# In the empirical csv file, an f value of 0 indicates a linearly infeasible system; an f value of 1 indicates a system that is 
# linearly feasible, but not integrally feasible; and an f value of 2 indicates an integrally feasible system.
#
# A full profile of an input program would be completed by running this script with every combination of modify and feasibility. 
# In some cases, these runs will each individually require a large amount of time, and/or a large amount of memory. This script 
# was originally designed around being run in a screen session and will regularly print status information to the command line if
# run as a separate process with the ampersand.
#
# Any program profiled with this script must print feasibility and timing information to stdout or stderr, with commas delimiting 
# separate values
#
# This script requires that all of these input files exist within an "input/" directory, and that "output/" and "empirical/" 
# directories also exist within the current working directory.

program="${1}"
if [[ ! -f "${program}" ]] ; then
  echo "Specified program, ${program}, not found."
  exit 1
fi
programNoExt="${program%%.*}"

modify="${2}"
if [[ "${modify}" = "variables" ]] ; then
  variablesArray=({1000..20000..1000})
  constraintsArray=(100000)
elif [[ "${modify}" = "constraints" ]] ; then
  variablesArray=(1000)
  constraintsArray=({10000..200000..10000})
else
  echo "Element to be modified must be either \"variables\" or \"constraints\"."
  exit 1
fi

f="${3}"
if [[ "${f}" != "f0" && "${f}" != "f1" && "${f}" != "f2" ]] ; then
  echo "Feasibility selection must be \"f0\", \"f1\", or \"f2\"."
  exit 1
fi

date=`date +%Y-%m-%d`

csvFile="empirical/${date}_${programNoExt}_Stuckey_${modify}_${f}.csv"

# For use with non-Opt* implementations:
echo "input file,variables,constraints,f,setup,linear,integer,cleanup,total,user,system,maximum resident set size (kB)" > "${csvFile}"

# For use with Opt* implementations:
#echo "input file,variables,constraints,f,false positives,main loop iterations,negative cycle length,setup,linear,integer,cleanup,total,user,system,maximum resident set size (kB)" > "${csvFile}"

for variables in "${variablesArray[@]}" ; do
  for constraints in "${constraintsArray[@]}" ; do
    for number in {1..3} ; do
      inputFile="input/test_n${variables}_m${constraints}_${f}_no${number}.utv"
      if [[ -f "${inputFile}" ]] ; then
        echo "${programNoExt} ${modify} ${f} processing ${inputFile}"
        outputFile="output/${date}_${programNoExt}_test_n${variables}_m${constraints}_${f}_no${number}_out.txt"
        echo -n "${inputFile},${variables},${constraints}," >> "${csvFile}"
        /usr/bin/time -f "%U,%S,%M" ./${program} "${inputFile}" "${outputFile}" &>> "${csvFile}" 
      fi
    done
  done
done
echo "${programNoExt} ${modify} ${f} complete"
