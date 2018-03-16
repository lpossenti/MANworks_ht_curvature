GETFEM_PREFIX=$(mkGetfemInc)/../
ifeq ($(WITH_SAMG),1)
CXXFLAGS += -I${SAMG}/
LFLAGS += -L/opt/lib/samg/
CXXFLAGS+= -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE -DPYRAMID_TRIANGULAR_FACETS
LIBS += -lamg -liomp5
endif

# getfem
CXXFLAGS+=$(shell getfem-config --cflags)
LDFLAGS+=$(shell getfem-config --libs)  

# superlu
# CXXFLAGS+=-DGMM_USES_SUPERLU -I$(mkSuperluInc)
# LDFLAGS+=-L$(mkSuperluLib) -lsuperlu
LDFLAGS+=-L$(mkQhullLib)

CXXFLAGS+=-std=c++14
CXX=g++-5
