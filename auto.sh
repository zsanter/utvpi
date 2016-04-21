#!/bin/bash

#./auto.sh [program] [modify:"variables"|"constraints"] [feasibility:"f0"|"f1"|"f2"]

program="${1}"
if [[ ! -f "${program}" ]] ; then
  echo "Specified program, ${program}, not found."
  exit 1
fi
modify="${2}"
if [[ "${modify}" != "variables" && "${modify}" != "constraints" ]] ; then
  echo "Element to be modified must be either \"variables\" or \"constraints\"."
  exit 1
fi
f="${3}"
if [[ "${f}" != "f0" && "${f}" != "f1" && "${f}" != "f2" ]] ; then
  echo "Feasibility selection must be \"f0\", \"f1\", or \"f2\"."
  exit 1
fi

DATE=`date +%Y-%m-%d`

csvFile="empirical/${program}_Stuckey_${f}_${modify}_${DATE}.csv"
echo "input file,variables,constraints,setup,linear,integer,cleanup,total,user,system,maximum resident set size (kb)" > "${csvFile}"

if [[ "${modify}" = "variables" ]] ; then
  constraints=100000
  for variables in {1000..20000..1000} ; do
    for number in {1..3} ; do
      inputFile="input/test_n${variables}_m${constraints}_${f}_no${number}.utv"
      if [[ -f "${inputFile}" ]] ; then
        echo "Processing ${inputFile}"
        outputFile="output/test_n${variables}_m${constraints}_${f}_no${number}_${DATE}_out.txt"
        echo -n "${inputFile},${variables},${constraints}," >> "${csvFile}"
        (/usr/bin/time -f "%U,%S,%M" ./${program} "${inputFile}" "${outputFile}") &>> "${csvFile}" 
      fi
    done
  done
elif [[ "${modify}" = "constraints" ]] ; then
  variables=1000
  for constraints in {10000..200000..10000} ; do
    for number in {1..3} ; do
      inputFile="input/test_n${variables}_m${constraints}_${f}_no${number}.utv"
      if [[ -f "${inputFile}" ]] ; then
        echo "Processing ${inputFile}"
        outputFile="output/test_n${variables}_m${constraints}_${f}_no${number}_${DATE}_out.txt"
        echo -n "${inputFile},${variables},${constraints}," >> "${csvFile}"
        (/usr/bin/time -f "%U,%S,%M" ./${program} "${inputFile}" "${outputFile}") &>> "${csvFile}" 
      fi
    done
  done
fi  