APP := integMask
CXX := g++
SRCDIR = ./src/
INCLDIR = ./includes/
BUILDDIR = ./build/
CXXFLAGS := -O3 -Wall -std=c++17
LDLIBS := -lstdc++fs -fopenmp -lhdf5_cpp -lhdf5 -lhdf5_hl -lhdf5_hl_cpp -lz -lcxsLib
SRC := $(wildcard $(SRCDIR)*.cpp)
OBJ := $(patsubst $(SRCDIR)%.cpp,$(BUILDDIR)%.o,$(SRC))
DEP := $(OBJ:.o=.d)

$(APP): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

-include $(DEP)

$(BUILDDIR)%.o: $(SRCDIR)%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INCLDIR) $(LDLIBS)
$(BUILDDIR)%.d: $(SRCDIR)%.cpp
	@$(CXX) $(CXXFLAGS) -I$(INCLDIR) $< -MM -MT $(@:.d=.o) >$@ $(LDLIBS)

init:
	mkdir -p $(BUILDDIR)

cleanAll: clean
	rm -f $(APP)

clean: 
	rm -f $(BUILDDIR)*.o $(BUILDDIR)*.d



