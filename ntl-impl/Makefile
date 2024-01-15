.PHONY: all

all: single parallel

single: main.cpp matvecop.hpp
	g++ -O3 main.cpp -DLOG_NAME='"log_single.csv"' -o single -lntl
parallel: main.cpp matvecop.hpp
	g++ -O3 main.cpp -DLOG_NAME='"log_parallel.csv"' -DDO_PARALLEL -o parallel -lntl
