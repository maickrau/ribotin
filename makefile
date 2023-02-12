GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Izstr/src -Iparallel-hashmap/parallel_hashmap/ -Wno-unused-parameter `pkg-config --cflags zlib`

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=`pkg-config --libs zlib`

_DEPS = fastqloader.h FastHasher.h VerkkoReadAssignment.h ReadExtractor.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = fastqloader.o FastHasher.o VerkkoReadAssignment.o ReadExtractor.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/seqpicker: $(OBJ) $(ODIR)/seqpicker.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/VerkkoIntegration: $(OBJ) $(ODIR)/VerkkoIntegration.o
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/seqpicker.o: $(SRCDIR)/seqpicker.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(ODIR)/VerkkoIntegration.o: $(SRCDIR)/VerkkoIntegration.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

all: $(BINDIR)/seqpicker $(BINDIR)/VerkkoIntegration

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
