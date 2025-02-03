GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Izstr/src -Iedlib/edlib/include -Iparallel-hashmap/parallel_hashmap/ -Icxxopts/include -Wno-unused-parameter `pkg-config --cflags zlib`
# silly workaround: bamtools does not have pkg-config cflags for finding the include directory
# instead assume it's a folder at the same location as zlib
CPPFLAGS+=`pkg-config --cflags zlib`/bamtools

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=`pkg-config --libs zlib` -lbamtools

_DEPS = fastqloader.h FastHasher.h VerkkoReadAssignment.h ReadExtractor.h KmerMatcher.h ClusterHandler.h VerkkoTangleGuesser.h RibotinUtils.h WfaHelper.h TwobitString.h TangleGuesser.h HifiasmIntegration.h CommonParams.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = fastqloader.o FastHasher.o VerkkoReadAssignment.o ReadExtractor.o KmerMatcher.o ClusterHandler.o VerkkoTangleGuesser.o RibotinUtils.o WfaHelper.o TwobitString.o TangleGuesser.o HifiasmIntegration.o CommonParams.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) $(LIBS) -lpthread -pthread -static-libstdc++

TEMPLATEPATH := $(realpath template_seqs)
VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

all: $(BINDIR)/ribotin-verkko $(BINDIR)/ribotin-hifiasm $(BINDIR)/ribotin-ref

$(BINDIR)/ribotin-verkko: $(OBJ) $(ODIR)/ribotin-verkko.o edlib/edlib/src/edlib.cpp
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ribotin-hifiasm: $(OBJ) $(ODIR)/ribotin-hifiasm.o edlib/edlib/src/edlib.cpp
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/ribotin-ref: $(OBJ) $(ODIR)/ribotin-ref.o edlib/edlib/src/edlib.cpp
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/CommonParams.o: $(SRCDIR)/CommonParams.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\"" -DRIBOTIN_TEMPLATE_PATH="\"$(TEMPLATEPATH)\""

$(ODIR)/ribotin-verkko.o: $(SRCDIR)/ribotin-verkko.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\"" -DRIBOTIN_TEMPLATE_PATH="\"$(TEMPLATEPATH)\""

$(ODIR)/ribotin-hifiasm.o: $(SRCDIR)/ribotin-hifiasm.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\"" -DRIBOTIN_TEMPLATE_PATH="\"$(TEMPLATEPATH)\""

$(ODIR)/ribotin-ref.o: $(SRCDIR)/ribotin-ref.cpp $(DEPS) $(OBJ)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\"" -DRIBOTIN_TEMPLATE_PATH="\"$(TEMPLATEPATH)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
