# ==============================================================================
# Paths (Update these to match your system installation paths)
# ==============================================================================
#Dogora
#AMREX_HOME ?= /home/mahesh/LBNL_Projects/DanMartinLDRDFY27/amrex/install
#HDF5_HOME  ?= /home/mahesh/Libraries/HDF5

#Perlmutter
AMREX_HOME ?= /pscratch/sd/n/nataraj2/DanMartinLDRDFY27/amrex/install
HDF5_HOME  ?= /opt/cray/pe/hdf5-parallel/1.14.3.1/gnu/12.3

# ==============================================================================
# Compiler & Flags
# ==============================================================================
# Use mpic++ or C++ compiler matching your AMReX build
CXX      = mpicxx
CXXFLAGS ?= -O3 -std=c++20 -Wall

# ==============================================================================
# Includes and Libraries
# ==============================================================================
# Add include directories for AMReX and HDF5
INCLUDES := -I. \
            -I$(AMREX_HOME)/include \
            -I$(HDF5_HOME)/include

# Add library search paths
LDFLAGS  := -L$(AMREX_HOME)/lib \
            -L$(HDF5_HOME)/lib

# Specify libraries to link against
# Note: Order matters. Dependent libraries come after.
LDLIBS   := -lamrex -lamrex_3d -lhdf5_hl -lhdf5 -lz -ldl -lm

# ==============================================================================
# Source and Target Definitions
# ==============================================================================
TARGET   := main.exe
SRCS     := main.cpp ReadChomboHDF5.cpp
OBJS     := $(SRCS:.cpp=.o)
DEPS     := ReadChomboHDF5.H

# ==============================================================================
# Build Rules
# ==============================================================================
.PHONY: all clean

all: $(TARGET)

# Link executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

# Compile object files
%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
