CXX=clang++
CXXFLAGS=-Ofast -g -Wall `root-config --cflags` -I$(ROOT_UTILS_PATH)/include -I$(PYTHIA8_DIR)/include -I$(ATLAS_PATH)/include
LDFLAGS=-Wl,-rpath,`root-config --libdir` `root-config --glibs` -L$(ROOT_UTILS_PATH)/lib -L$(PYTHIA8_DIR)/lib -L$(ATLAS_PATH)/lib -lUtilities -lAtlasStyle -ldl

reqdirs= bin output log errors Plots

directories:
	mkdir -p ${reqdirs}

all: directories gen analyze

gen: src/gen.cxx
	$(CXX) $< $(CXXFLAGS) $(LDFLAGS) $(PYTHIA8_DIR)/lib/libpythia8.a -o bin/gen

analyze: src/analyze.cxx
	$(CXX) $< $(CXXFLAGS) $(LDFLAGS) -o bin/analyze

clean:
	rm -rf bin logs errors
