APP := ../libtinyxml.a
CXX := g++
SRCDIR = ./src/
INCLDIR = ./includes/
BUILDDIR = ./build/
CXXFLAGS := -O3 -Wall -std=c++17 
SRC := $(wildcard $(SRCDIR)*.cpp)
OBJ := $(patsubst $(SRCDIR)%.cpp,$(BUILDDIR)%.o,$(SRC))
DEP := $(OBJ:.o=.d)

$(APP): $(OBJ)
	ar rcs $@ $^

-include $(DEP)

$(BUILDDIR)%.o: $(SRCDIR)%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INCLDIR) 
$(BUILDDIR)%.d: $(SRCDIR)%.cpp
	@$(CXX) $(CXXFLAGS) -I$(INCLDIR) $< -MM -MT $(@:.d=.o) >$@ 

init:
	mkdir -p $(BUILDDIR)

cleanAll: clean
	rm -f $(APP)

clean: 
	rm -f $(BUILDDIR)*.o $(BUILDDIR)*.d



