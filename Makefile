MAKEFILE      = Makefile

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -std=c++14 # -Wextra -pedantic
CFLAGS        = -O2 -Wall $(DEFINES)
CXXFLAGS      = -O2 -Wall $(DEFINES)
INCPATH       = -I.
LINK          = g++
LFLAGS        = -Wl,-O1
LIBS          = $(SUBLIBS) # -lpthread -lstdc++
AR            = ar cqs
RANLIB        = 
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = 
OBJECTS       = support.o	mathmore.o	main.o	main_test.o	commandline.o	## TODO: add all corresponding .o files here
DIST          =
DESTDIR       = #avoid trailing-slash linebreak
TARGET        = support		test


first: all
####### Implicit rules

# .SUFFIXES: .o .c .cpp .cc .cxx .C

# .cpp.o:
# 	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

# .cc.o:
# 	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

# .cxx.o:
# 	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

# .C.o:
# 	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

# .c.o:
# 	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: $(TARGET)

support: support.o main.o commandline.o mathmore.o
	$(LINK) $(LFLAGS) support.o main.o commandline.o mathmore.o $(LIBS) -o support

test: support.o main_test.o mathmore.o
	$(LINK) $(LFLAGS) support.o mathmore.o main_test.o $(LIBS) -o test

# $(TARGET):  $(OBJECTS)  
# 	$(LINK) $(LFLAGS) $(OBJECTS) $(OBJCOMP) $(LIBS) -o $(TARGET)

dist: 
	@test -d .tmp/test1.0.0 || mkdir -p .tmp/test1.0.0
	$(COPY_FILE) --parents $(DIST) .tmp/test1.0.0/ && $(COPY_FILE) --parents hello.h .tmp/test1.0.0/ && $(COPY_FILE) --parents hello.cpp main.cpp .tmp/test1.0.0/ && (cd `dirname .tmp/test1.0.0` && $(TAR) test1.0.0.tar test1.0.0 && $(COMPRESS) test1.0.0.tar) && $(MOVE) `dirname .tmp/test1.0.0`/test1.0.0.tar.gz . && $(DEL_FILE) -r .tmp/test1.0.0


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)

distclean: clean 
	-$(DEL_FILE) $(TARGET) 

####### Sub-libraries

check: first

compiler_clean: 

####### Compile

# TODO: add all dependencies for the .o

commandline.o: commandline.cpp commandline.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o commandline.o commandline.cpp

support.o: support.cpp support.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o support.o support.cpp

main.o: main.cpp support.h commandline.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

mathmore.o: mathmore.cpp mathmore.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o mathmore.o mathmore.cpp

main_test.o: main_test.cpp support.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main_test.o main_test.cpp
