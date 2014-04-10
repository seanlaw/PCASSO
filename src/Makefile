#Sean M. Law

include ../makefile.defs

OBJS = libpcasso.$(LIBEXT)

all:	$(BINDIR) libpcasso.$(LIBEXT) pcasso 

clean:
	rm -rf *.o
	rm -rf *stackdump
	rm -rf ../bin/*
	rm -rf *~

$(BINDIR):
	mkdir -p $(BINDIR)

libpcasso.$(LIBEXT):
	$(MAKE) -C ../lib libpcasso.$(LIBEXT)

pcasso: pcasso.o $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) $(LINKFLAGS) -o $(BINDIR)/$@ $<

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -c $(INCFLAGS) $<
