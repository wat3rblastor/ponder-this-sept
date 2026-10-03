CC      ?= cc
CFLAGS  ?= -O3 -march=native -funroll-loops -Wall -Wextra -std=gnu11
BIN     := build/loeschsearch

.PHONY: check test bin clean

check: test bin
	@echo "=== make check: OK ==="

test:
	python3 -m unittest discover -s tests -t . -v

bin: $(BIN)

$(BIN): src/c/loeschsearch.c | build
	$(CC) $(CFLAGS) -o $@ $<

build:
	mkdir -p build

clean:
	rm -rf build

# CUDA engine (NVIDIA). Requires nvcc; sm arch from the installed GPU.
NVCC     ?= nvcc
CUDA_ARCH ?= sm_121
cuda: build/apsearch_cuda
build/loesch_core_api.o: src/c/loesch_core_api.c src/c/loesch_core.h | build
	gcc -O3 -march=native -std=gnu11 -c -o $@ $<
build/apsearch_cuda: src/c/apsearch_cuda.cu src/c/loesch_core_api.h build/loesch_core_api.o | build
	$(NVCC) -O3 -arch=$(CUDA_ARCH) -Xcompiler -fopenmp -I src/c -o $@ src/c/apsearch_cuda.cu build/loesch_core_api.o -lm
apsearch: build/apsearch
build/apsearch: src/c/apsearch.c src/c/loesch_core.h | build
	$(CC) -O3 -march=native -std=gnu11 -o $@ $< -lm
