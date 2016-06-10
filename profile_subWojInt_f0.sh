#!/bin/bash

./auto.sh subWojInt.exe variables f0 &
./auto.sh subWojInt.exe constraints f0 &

wait
echo "f0 profiling complete."
