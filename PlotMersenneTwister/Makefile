CC=g++
CFLAGS=-Wall `root-config --cflags`
LIBS=`root-config --libs`

PlotMersenneTwister: PlotMersenneTwister.o
	$(CC) $(LIBS) -o $@ $<

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@
