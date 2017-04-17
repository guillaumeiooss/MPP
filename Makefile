OBJS = linAlg.o MPP_rect.o test_MPP.o
BIN = MPP

all: $(OBJS)
	g++ $(OBJS) -Wall -o $(BIN)

%.o: %.cpp
	g++ -Wall -c $<

clean:
	-rm -f $(OBJS) *~

realclean:
	$(MAKE) clean
	-rm -f $(BIN)

