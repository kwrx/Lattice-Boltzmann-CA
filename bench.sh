#!/bin/sh


run_wind() {


    export NPROCS=1


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=100" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r1=$(cat /tmp/bench-output)


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=1000" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r2=$(cat /tmp/bench-output)


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=10000" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r3=$(cat /tmp/bench-output)




    export NPROCS=2


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=100" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r4=$(cat /tmp/bench-output)


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=1000" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r5=$(cat /tmp/bench-output)


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=10000" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r6=$(cat /tmp/bench-output)




    export NPROCS=4


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=100" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r7=$(cat /tmp/bench-output)


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=1000" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r8=$(cat /tmp/bench-output)


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=10000" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r9=$(cat /tmp/bench-output)




    export NPROCS=8


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=100" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r10=$(cat /tmp/bench-output)


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=1000" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r11=$(cat /tmp/bench-output)


    make -s clean
    make -s CXXFLAGS="-DBENCH -DITERATIONS=10000" run

    mpirun -np $NPROCS --oversubscribe $1 &> /tmp/bench-output
    r12=$(cat /tmp/bench-output)


    ./graph.py 100 1000 10000 $r1 $r2 $r3 $r4 $r5 $r6 $r7 $r8 $r9 $r10 $r11 $r12

}


run_wind apsd