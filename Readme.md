# Artifacts for "Faster Signatures from MPC-in-the-Head"

This folder contains the artifacts for the paper "Faster Signatures from MPC-in-the-Head" accepted at ASIACRYPT 2024. This includes the necessary source code to reproduce the results from the paper.

The artifact contains two folders: 
    - ```ref```: a reference implementation of the signature specifically tailored to be efficient in the presence of large memory 
    - ```poc```: an implementation made to be efficient using less memory but also containing the different ways to compare the timing with the SDitH project.

The code was primarily developed on Intel/AMD processors with AES instructions available. We currently do not support other platforms.

To use the code it is necessary to have OpenSSL and BSD message digest installed as the development package, which can be done as follows for instance:
```sudo apt-get install build-essential libmd-dev``` (for ```ref```)
```sudo apt-get install libssl-dev``` (for ```poc```)



## Reference implementation 

The code from 'ref' contains the implementations of three versions of the signature: a version with a tightly optimized folding, a version with classical folding, and a version with classical folding and the use of SHA instead of AES.
All versions use the 128-bit security level with parameter set K = 1736, k = 960, w = 217, bs = 8 as specified in Section 7.

To use the more portable version, it is necessary to install crypto-algorithms by Brad Conte. This is achieved by running the command:
```git clone https://github.com/B-Con/crypto-algorithms.git```
in the root directory of the present code artifact.

The number of parties used in the signature scheme can be specified using the '-DL$i' command. We use the parameters for the 128-bit security level specified in Table 1 on page 18 of the paper. The number of parties and the number of rounds is set according the value of $i given.
$i corresponds to the log of the number of parties N=2^D that should be used.
The different versions using '-DL8', '-DL9', '-DL10', '-DL11', '-DL12', '-DL13', '-DL15', '-DL16' correspond to the parameters for D and tau specified in Table 3 in the paper.

The option '-DSANITIZE' allows to set the memory of sensitive data to zero when they are not used anymore.

The option '-DOPTIM_SCALARS' allows using fewer scalar multiplications by making full use of the register size.

The option '-DCOMPRESSION' version that uses an optimized folding for the GGM tree which requires extra work memory size.

The option '-DCLASSICAL_SHA_TREE' uses SHA256 in the construction of the GGM tree. If this is not specified, a construction based on AES is used.

The option '-DCLASSICAL_FOLDING' uses the standard folding of the tree. If this is not specified, the fast folding algorithm described in Section 6.2 is used in the protocol.

By running the file 'maketest.sh', a script is run which executes every parameter choice for the different versions of the signature specified above (fast/classical folding, GGM tree with SHA).
The outputs include the work memory size, the signature size, as well as the timings and clock cycles of the key generation, signature and verification algorithms.

To run the portable version, execute 'maketest.sh portable' to use the portable AES and SHA functions from the "crypto-algorithms" repository.

Warning:
The Intel/AMD versions and the portable versions do not give the exact same keys or signatures, due to a different endianness for the respective AES implementations.


## POC implementation

Proof-of-concept implementation of the signature using prior version of the parameters (namely it is specifically tailored for 16-bit blocksize). It is slower than the reference implementation but requires significantly less memory. It still offers better performances than sha-based GGM tree, for which the code also provides a comparison with the SDitH project. It is also the base used to be injected into the SditH project in order to have proper comparison of the impact of the AES-based tree in this project.

Several pieces of code as well as the ```sdith_hypercube_cat1_gf256``` folder have been extracted from the submission package of SDitH in order to have the exact same ways of computing benchmark and have accurate comparison between the two signature scheme.

The signature can be compile using the Makefile by simply running the command ```make```. By running ```bench``` with the number of signature to be generated and verified (e.g ```./bench 100``` for 100 tests), the signature config will be output on top of the number of correct verification (namely Verif(Sign(m))== True), the timing in milliseconds and clock cycles, and the communication costs. 

The modified files in the SDitH submission package in order to have a GGM tree that also support AES are the following:

    * rng.h: define some flag and the AES encryption functions.
    * sdith.c: modification to SDitH GGM tree construction to use AES instead.
    * param.h: flag to be present for AES to be unused. Uncomment to switch to AES mode.
    * treeprg.c: AES-based seed expand.