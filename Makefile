CC      ?= cc
CFLAGS  ?= -O3 -march=native -funroll-loops -Wall -Wextra -std=c11
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
