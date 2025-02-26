#!/bin/sh

if [ $# -ne 0 ]
then
    if [ "$1" = "portable" ]
    then
    echo "Running the portable version"
    echo "============================"
    DEFS="-DCOMPRESSION -DSANITIZE -DGENERICAES -DGENERICSHA"
    SRC="RSDhypercube.c RSDtools.c ../crypto-algorithms/aes.c ../crypto-algorithms/sha256.c"
    LIBINCLUDE="-lm"
    #loop="10 9 8"
    loop="8 9 10"
    fi
    if [ "$1" = "portablesha" ]
    then
    echo "Running the Intel/AMD version without libmd"
    echo "============================"
    DEFS="-DCOMPRESSION -DSANITIZE -DGENERICSHA"
    SRC="RSDhypercube.c RSDtools.c ../crypto-algorithms/sha256.c"
    LIBINCLUDE="-lm"
    fi
else
    echo "Running the Intel/AMD version"
    echo "============================="
    DEFS="-DCOMPRESSION -DOPTIM_SCALARS -DSANITIZE"
    SRC="RSDhypercube.c RSDtools.c"
    LIBINCLUDE="-lm -lmd"
fi

loop="16 15 13 12 11 10 9 8"

echo "FAST FOLDING"
for i in $loop; do echo $i;gcc -O3 -DL$i $DEFS -march=native -o RSDhypercube $SRC $LIBINCLUDE;  ./RSDhypercube ; done
echo "----------------------------"

echo "SLOW FOLDING"
for i in $loop; do echo $i;gcc -O3 -DL$i -DCLASSICAL_FOLDING $DEFS -march=native -o RSDhypercube $SRC $LIBINCLUDE;  ./RSDhypercube ; done
echo "----------------------------"

echo "SHA256 and SLOW FOLDING"
for i in $loop; do echo $i;gcc -O3 -DL$i -DCLASSICAL_SHA_TREE -DCLASSICAL_FOLDING $DEFS -march=native -o RSDhypercube $SRC $LIBINCLUDE;  ./RSDhypercube ; done
echo "----------------------------"
