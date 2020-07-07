
OBJSLIB = linAlg.o MPP_rect.o lexmin.o MPP_gen.o
OBJS = $(OBJSLIB) test_MPP.o
BIN = MPP.out
LIB = libmpp.a
LIBSO = libmpp.so

all: $(OBJS)
	g++ $(OBJS) -Wall -L/usr/local/lib -lpiplib64 -o $(BIN)

%.o: %.cpp
	g++ -Wall -c -fPIC $<

lib: $(OBJSLIB)
	ar rcs $(LIB) $^

static: $(OBJS)
	g++ $(OBJS) -shared -Wall -L/usr/local/lib -lpiplib64 -o $(LIBSO)

clean:
	-rm -f $(OBJS) *~

realclean:
	$(MAKE) clean
	-rm -f $(BIN) $(LIB)

