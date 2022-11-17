#!/bin/bash
cd ..
make clean
make devel
cd ./tests
../bin/mdlovofit -f 0.7 -t output.pdb ./data/trajectory.pdb > rmsd.dat
if [ ! -z `diff -I"Version" rmsd.dat ./results/TEST1.dat` ]; then
    echo "ERROR in TEST1"
    exit
else
    echo "TEST1 passed!"
fi
../bin/mdlovofit -iref 13 -f 0.7 -t output.pdb ./data/trajectory.pdb > rmsd.dat
if [ ! -z `diff -I"Version" rmsd.dat ./results/TEST2.dat` ]; then
    echo "ERROR in TEST2"
    exit
else
    echo "TEST2 passed!"
fi
split --lines 8784 ./data/trajectory.pdb split
../bin/mdlovofit -iref 13 -f 0.7 -t output.pdb splitaa splitab splitac splitad > rmsd.dat
if [ ! -z `diff -I"Version" rmsd.dat ./results/TEST3.dat` ]; then
    echo "ERROR in TEST3"
    exit
else
    echo "TEST3 passed!"
fi
\rm -f rmsd.dat output.pdb split*
