#!/bin/bash

program="${1}"
if [[ ! -f "${program}" ]] ; then
  echo "Specified program, ${program}, not found."
  exit 1
fi

"${program}" "variables" "f2"
"${program}" "constraints" "f2"
"${program}" "variables" "f1"
"${program}" "constraints" "f1"
"${program}" "variables" "f0"
"${program}" "constraints" "f0"
