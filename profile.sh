#!/bin/bash

program="${1}"
if [[ ! -f "${program}" ]] ; then
  echo "Specified program, ${program}, not found."
  exit 1
fi

./auto.sh "${program}" "variables" "f2"
./auto.sh "${program}" "constraints" "f2"
./auto.sh "${program}" "variables" "f1"
./auto.sh "${program}" "constraints" "f1"
./auto.sh "${program}" "variables" "f0"
./auto.sh "${program}" "constraints" "f0"
