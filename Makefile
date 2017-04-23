OBJS = linAlg.o MPP_rect.o lexmin.o test_MPP.o
BIN = MPP.out
# MPP_gen.o

all: $(OBJS)
	g++ $(OBJS) -Wall -L/usr/local/lib -lpiplib64 -o $(BIN)

%.o: %.cpp
	g++ -Wall -c $<

clean:
	-rm -f $(OBJS) *~

realclean:
	$(MAKE) clean
	-rm -f $(BIN)

