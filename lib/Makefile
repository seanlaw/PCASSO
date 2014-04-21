#Sean M. Law

include ../makefile.defs

LIBOBJS = Analyze.o Atom.o Chain.o DTree.o Misc.o Molecule.o PCASSO.o PDB.o Prmtop.o Residue.o Select.o Trajectory.o Coor.o 

all:	libpcasso.$(LIBEXT)

libpcasso.$(LIBEXT):	$(LIBOBJS)
	ar rc libpcasso.$(LIBEXT) $(LIBOBJS)

clean:
	rm -rf *.o
	rm -rf *.$(LIBEXT)
	rm -rf *stackdump
	rm -rf *~

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -c $< -o $@
