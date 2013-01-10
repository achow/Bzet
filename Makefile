GCC = g++
CFLAGS =-g -Wall -DBZET_IMPL_ -DNTESTS=10 #-DDEBUG
LIBS =

four: Bzet.h tester.cpp
	$(GCC) $(CFLAGS) $(LIBS) -o Bzet4Test -DNODE_ELS=4 Bzet.h tester.cpp
	./Bzet4Test

eight: Bzet.h tester.cpp
	$(GCC) $(CFLAGS) $(LIBS) -o Bzet8Test -DNODE_ELS=8 Bzet.h tester.cpp
	./Bzet8Test

all: Bzet.h tester.cpp
	$(GCC) $(CFLAGS) $(LIBS) -o Bzet4Test -DNODE_ELS=4 Bzet.h tester.cpp
	$(GCC) $(CFLAGS) $(LIBS) -o Bzet8Test -DNODE_ELS=8 Bzet.h tester.cpp
	$(GCC) $(CFLAGS) $(LIBS) -o Bzet16Test -DNODE_ELS=16 Bzet.h tester.cpp
	$(GCC) $(CFLAGS) $(LIBS) -o Bzet32Test -DNODE_ELS=32 Bzet.h tester.cpp
	$(GCC) $(CFLAGS) $(LIBS) -o Bzet64Test -DNODE_ELS=64 Bzet.h tester.cpp

check: all
	./Bzet4Test && ./Bzet8Test && ./Bzet16Test && ./Bzet32Test && ./Bzet64Test

dummy:
	./Bzet4Test || echo 'Bzet4 test failed'
	./Bzet8Test || echo 'Bzet8 test failed'
	./Bzet16Test || echo 'Bzet16 test failed'
	./Bzet32Test || echo 'Bzet32 test failed' 
	./Bzet64Test || echo 'Bzet64 test failed'
