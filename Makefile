PROG=mlinker
UNAME_S := $(shell uname -s)

# CC
ifeq ($(UNAME_S),Darwin)
  CC := clang++ -std=c++11 -stdlib=libc++
else
  CC := c++ -o3 -std=c++11
endif

LDIR=./packages/bamtools/lib/
IDIR2=./htslib/
LDIR2=./packages/htslib/lib/
IDIR=./packages/bamtools/include/bamtools/
IDIR3=./packages/bamtools/lib/

CFLAGS=-I$(IDIR) -I$(IDIR2)
LFLAGS=-L$(LDIR) -L$(LDIR2) -Wl,-rpath,$(IDIR2),-rpath,$(IDIR3) #-L$(LDIR3)
LBAMTOOLS=-lbamtools
LHTS=-lhts
LZLIB=-lz
LCURL=-lcurl
LIBRARIES=$(LHTS) $(LBAMTOOLS) $(LCURL) $(LZLIB)


#CCPLUS=$(CC) $(LFLAGS) $(CFLAGS)

SRCDIR:=src
#BINDIR:=bin
#TARGET:=./$(BINDIR)/$(PROG)
TARGET:=./$(PROG)
SOURCES=$(shell find ./$(SRCDIR)/*.cpp)
OBJECTS=$(addsuffix .o,$(basename $(SOURCES)) )

##########################
all: $(TARGET)

##########################
$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(TARGET) $(LFLAGS) $(CFLAGS) $(LIBRARIES)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(LFLAGS) $(CFLAGS) -c -o $@ $<

print-% : ; @echo $* = $($*)

##########################
clean:
	-@rm $(OBJECTS)
	-@rm -v $(TARGET)

distclean:
	-@rm $(OBJECTS)
	-@rm -v $(TARGET)
	-@rm ./output/variant_network_*
	-@rm ./output/hap_solution_*
	-@rm ./output/het_coverage_*
	-@rm ./output/bxmatrix_*
	-@rm ./output/filtered_coverage_*
	-@rm ./output/sv_phased_*
	-@rm ./output/cn_phased_*
	-@rm ./output/hic_links_*
	-@rm ./output/phased_sv_sites.dat

##########################
batch:
	make clean
	make all
