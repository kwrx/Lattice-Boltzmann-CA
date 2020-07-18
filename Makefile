.PHONY: all clean


OUTPUT 	:= apsd
SRCS	:= src/main.cpp

OPT 	:= -O3 -g -fno-stack-protector
LIBS	:= -lallegro -lallegro_primitives -lallegro_font


CXX		:= mpic++
NPROCS	?= 2


all: $(OUTPUT)
$(OUTPUT): $(SRCS)
	$(CXX) $(CXXFLAGS) $(OPT) -o $@ $< $(LIBS)

bench:
	./bench.sh
	
clean:
	$(RM) $(OUTPUT)

debug: $(OUTPUT)
	chmod +x $<
	mpirun -np $(NPROCS) --oversubscribe --display-map xterm -e gdb -ex run $<

run: $(OUTPUT)
	chmod +x $<
	mpirun -np $(NPROCS) --oversubscribe $<


run-nompi: $(OUTPUT)
	chmod +x $<
	./$<
