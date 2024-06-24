GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Izstr/src -Iedlib/edlib/include -Iparallel-hashmap/parallel_hashmap/ -Icxxopts/include -Wno-unused-parameter `pkg-config --cflags zlib` -Iconcurrentqueue -Ihifioverlapper/src -Ifastacompressor/src

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=`pkg-config --libs zlib`

_DEPS = Common.h UnionFind.h TwobitString.h SparseEdgeContainer.h RankBitvector.h VectorWithDirection.h MostlySparse2DHashmap.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = Common.o UnionFind.o TwobitString.o SparseEdgeContainer.o RankBitvector.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/chunkgraph: $(OBJ) $(ODIR)/chunkgraph.o edlib/edlib/src/edlib.cpp fastacompressor/lib/fastacompress.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

fastacompressor/lib/fastacompress.a:
	$(MAKE) -C fastacompressor lib

all: $(BINDIR)/chunkgraph

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
	$(MAKE) -C fastacompressor clean
