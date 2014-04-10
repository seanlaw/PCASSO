#Sean M. Law

include ../makefile.defs

LIBOBJS = Analyze.o Atom.o BTree.o Chain.o Constants.o DTree.o Enum.o Misc.o Molecule.o PCASSO.o PDB.o Prmtop.o Residue.o Select.o Trajectory.o Vector.o 

all:	libpcasso.$(LIBEXT)

libpcasso.$(LIBEXT):	$(LIBOBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $^ -o libpcasso.$(LIBEXT)

clean:
	rm -rf *.o
	rm -rf *.$(LIBEXT)
	rm -rf *stackdump
	rm -rf *~

%.o:	%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $(DEBUG) -c $< -o $@
