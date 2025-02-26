# This Makefile was inspired by SDITH submission package

.POSIX:

HASH_PATH=sha3
HASH_INCLUDE=-Isha3/avx2

ifndef $(HOST_OS)
	ifeq ($(OS),Windows_NT)
		HOST_OS = Windows
	else
		HOST_OS := $(shell uname -s)
	endif
endif

ifeq ($(HOST_OS),Darwin)
$(error Darwin platform not supported yet)
endif

AVXFLAGS = -DAVX2 -mavx2 -mpclmul -mgfni -mavx -maes -mpclmul -msse4.2 -mrdrnd
HASH_MAKE_OPTIONS = PLATFORM=avx2 

CC = gcc
CFLAGS = -W -Wall -O3 -fPIC -std=c11 -flto

CRYPTO_SRC = hash-sha3.c pprf_sig.c
BENCH_MAIN_SRC = benchmark/bench.c benchmark/timing.c

all: bench

libhash:
	$(HASH_MAKE_OPTIONS) make -C $(HASH_PATH)

sign: ${SRC} ${CRYPTO_SRC} libhash
	-${CC} ${CFLAGS} ${AVXFLAGS} -o $@ ${SRC} ${CRYPTO_SRC} -I. ${HASH_INCLUDE} -L${HASH_PATH} -lcrypto -lhash

bench: ${SRC} ${CRYPTO_SRC} ${BENCH_MAIN_SRC} libhash
	-${CC} ${CFLAGS} ${AVXFLAGS} -DBENCHMARK_CYCLES -o $@ ${SRC} ${CRYPTO_SRC} ${BENCH_MAIN_SRC} -I. ${HASH_INCLUDE} -L${HASH_PATH} -lcrypto -lhash -lm

clean:
	-rm -f *.o *.a sign bench
	$(HASH_MAKE_OPTIONS) make -C $(HASH_PATH) clean