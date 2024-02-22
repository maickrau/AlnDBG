GPP=$(CXX)
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -g -Izstr/src -Iparallel-hashmap/parallel_hashmap/ -Icxxopts/include -Wno-unused-parameter `pkg-config --cflags zlib` -IMBG/src -Iconcurrentqueue -Ihifioverlapper/src

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=`pkg-config --libs zlib`

_DEPS = KmerGraph.h KmerMatcher.h MatchGroup.h UnitigGraph.h GraphCleaner.h AlnHaploFilter.h GraphPhaser.h Common.h GraphResolver.h AnchorFinder.h UnionFind.h ChunkmerFilter.h MultiplexResolverCaller.h
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ = KmerGraph.o KmerMatcher.o UnitigGraph.o GraphCleaner.o AlnHaploFilter.o GraphPhaser.o Common.o GraphResolver.o AnchorFinder.o UnionFind.o ChunkmerFilter.o MultiplexResolverCaller.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/AlnDBG: $(OBJ) $(ODIR)/main.o hifioverlapper/lib/hifioverlapper.a MBG/lib/mbg.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/chunkgraph: $(OBJ) $(ODIR)/chunkgraph.o hifioverlapper/lib/hifioverlapper.a MBG/lib/mbg.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(BINDIR)/alncorrect: $(OBJ) $(ODIR)/alncorrect.o hifioverlapper/lib/hifioverlapper.a MBG/lib/mbg.a
	$(GPP) -o $@ $^ $(LINKFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -c -o $@ $< $(CPPFLAGS)

MBG/lib/mbg.a:
	$(MAKE) -C MBG lib

hifioverlapper/lib/hifioverlapper.a:
	$(MAKE) -C hifioverlapper lib

all: $(BINDIR)/AlnDBG $(BINDIR)/alncorrect $(BINDIR)/chunkgraph

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
	$(MAKE) -C MBG clean
	$(MAKE) -C hifioverlapper clean
