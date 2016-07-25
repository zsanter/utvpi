#!/bin/bash

#./auto.sh [program] [modify:"variables"|"constraints"] [feasibility:"f0"|"f1"|"f2"]

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
echo "input file,variables,constraints,f,setup,linear,integer,cleanup,total,user,system,maximum resident set size (kB)" > "${csvFile}"
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
