CXXC=g++
LIBS=-lm -lz
CFLAGS = -O3 -g
HG_DEFS = -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
HG_WARN=-Wformat -Wreturn-type
BIO_DIR = ./bioUtils
INCLUDES = -I$(BIO_DIR)
BIO_LIBS = -L$(BIO_DIR)/ -lbiotools

all:
	cd $(BIO_DIR); make
	make kturnSeeker

clean:
	cd $(BIO_DIR); make clean
	rm -f *.o

kturnSeeker: kturnSeeker.o kturnSeekerMain.o
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -o kturnSeeker kturnSeekerMain.o kturnSeeker.o $(BIO_LIBS) $(LIBS)

kturnSeeker.o: kturnSeeker.cpp kturnSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c kturnSeeker.cpp
  
kturnSeekerMain.o: kturnSeekerMain.cpp kturnSeeker.h
	$(CXXC) $(CFLAGS) ${HG_DEFS} ${HG_WARN} $(INCLUDES) -c kturnSeekerMain.cpp



