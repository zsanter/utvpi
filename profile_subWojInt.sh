#!/bin/bash

./auto.sh subWojInt.exe variables f1 &
./auto.sh subWojInt.exe variables f2 &
./auto.sh subWojInt.exe constraints f1 &
./auto.sh subWojInt.exe constraints f2 &

wait
echo "All profiling complete."
