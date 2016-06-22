#!/bin/bash

./auto.sh lahiri.exe variables f1 &
./auto.sh lahiri.exe variables f2 &
./auto.sh lahiri.exe constraints f1 &
./auto.sh lahiri.exe constraints f2 &

wait
echo "All profiling complete."
