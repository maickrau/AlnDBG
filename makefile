GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Izstr/src -Iedlib/edlib/include -Iparallel-hashmap/parallel_hashmap/ -Icxxopts/include -Wno-unused-parameter `pkg-config --cflags zlib` -Iconcurrentqueue -Ifastacompressor/src

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=`pkg-config --libs zlib`

_DEPS = Common.h UnionFind.h TwobitString.h SparseEdgeContainer.h RankBitvector.h ChunkGraphWriter.h ChunkUnitigGraph.h EdlibWrapper.h ConsensusMaker.h KmerIterator.h GraphCleaner.h SequenceHelper.h ChunkExtractor.h SequenceIdentitySplitter.h ChunkPhasing.h ChunkResolution.h CanonHelper.h OverlapMatcher.h PathWalker.h TransitiveClosure.h MinimizerCorrection.h ChunkHelper.h VectorWithDirection.h MostlySparse2DHashmap.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = Common.o UnionFind.o TwobitString.o SparseEdgeContainer.o RankBitvector.o ChunkGraphWriter.o ChunkUnitigGraph.o EdlibWrapper.o ConsensusMaker.o KmerIterator.o GraphCleaner.o SequenceHelper.o ChunkExtractor.o SequenceIdentitySplitter.o ChunkPhasing.o ChunkResolution.o CanonHelper.o OverlapMatcher.o PathWalker.o TransitiveClosure.o MinimizerCorrection.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/chunkgraph: $(OBJ) $(ODIR)/chunkgraph.o edlib/edlib/src/edlib.cpp fastacompressor/lib/fastacompress.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/chunkgraph.o: $(SRCDIR)/chunkgraph.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

fastacompressor/lib/fastacompress.a:
	$(MAKE) -C fastacompressor lib

all: $(BINDIR)/chunkgraph

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
	$(MAKE) -C fastacompressor clean
