BIN_DIR = bin/

all: ${BIN_DIR} dpc2sim-stream dpc2sim-ampm-lite build-dpc2

$(BIN_DIR):
	mkdir -p $@

run: dpc2sim-stream dpc2sim-ampm-lite
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-stream
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-ampa-lite

dpc2sim-stream:
	$(CXX) -Wall -no-pie -o bin/dpc2sim-stream example_prefetchers/stream_prefetcher.cc lib/dpc2sim.a

dpc2sim-ampm-lite:
	$(CXX) -Wall -no-pie -o bin/dpc2sim-ampa-lite example_prefetchers/ampm_lite_prefetcher.cc lib/dpc2sim.a

build-dpc2:
	$(CXX) -Wall -no-pie -o bin/dpc2sim-sbooe example_prefetchers/dpc2_brown.cpp lib/dpc2sim.a
	$(CXX) -Wall -no-pie -o bin/dpc2sim-oampm example_prefetchers/dpc2_jia.cpp lib/dpc2sim.a
	$(CXX) -Wall -no-pie -o bin/dpc2sim-spad  example_prefetchers/dpc2_karsli.c lib/dpc2sim.a
	$(CXX) -Wall -no-pie -o bin/dpc2sim-spp   example_prefetchers/dpc2_kim.cpp lib/dpc2sim.a
	$(CXX) -Wall -no-pie -o bin/dpc2sim-bop   example_prefetchers/dpc2_michaud.c lib/dpc2sim.a
	$(CXX) -Wall -no-pie -o bin/dpc2sim-vldp  example_prefetchers/dpc2_shevgoor.c lib/dpc2sim.a
	$(CXX) -Wall -no-pie -o bin/dpc2sim-tempo example_prefetchers/dpc2_sutherland.cc lib/dpc2sim.a
	$(CXX) -Wall -no-pie -o bin/dpc2sim-slim-ampm example_prefetchers/dpc2_young.c lib/dpc2sim.a

run-dpc2sim:
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-sbooe
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-oampm
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-spad
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-spp
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-bop
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-vldp
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-tempo
	zcat traces/mcf_trace2.dpc.gz | ./bin/dpc2sim-slim-ampm




clean:
	rm -rf bin/*

.PHONY: all run clean
