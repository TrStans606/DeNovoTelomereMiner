# Makefile adapted from http://www.partow.net/programming/makefile/index.html

CXX     := c++
CXXFLAGS    := -std=c++12 -pedantic-errors -Wall -Wextra #-Werror
LDFLAGS     := -L/usr/lib -lboost_system -lboost_program_options -lboost_filesystem -lboost_date_time -lstdc++ -lm
BUILD       := ./build
OBJ_DIR     := $(BUILD)/objects
SRC_DIR     := ./src
APP_DIR     := $(BUILD)/apps
INCLUDE     := -Iinclude/

TF_SRC      := $(shell find ./src/core ./src/telomereFinder -type f -name '*.cpp') ./src/telomereFinder_driver.cpp
TF_OBJ      := $(patsubst ./src/%, $(OBJ_DIR)/%, $(TF_SRC:.cpp=.o))
JF_SRC      := $(shell find ./src/core ./src/junctionFinder -type f -name '*.cpp') ./src/junctionFinder_driver.cpp
JF_OBJ      := $(patsubst ./src/%, $(OBJ_DIR)/%, $(JF_SRC:.cpp=.o))
SQ_SRC      := $(shell find ./src/core -type f -name '*.cpp') ./src/sequenceQuality_driver.cpp
SQ_OBJ      := $(patsubst ./src/%, $(OBJ_DIR)/%, $(SQ_SRC:.cpp=.o))
WCD_SRC     := $(shell find ./src/core -type f -name '*.cpp') ./src/wcdInterrogate_driver.cpp
WCD_OBJ     := $(patsubst ./src/%, $(OBJ_DIR)/%, $(WCD_SRC:.cpp=.o))

all: build telomereFinder junctionFinder sequenceQuality wcdInterrogate

telomereFinder: $(APP_DIR)/telomereFinder
junctionFinder: $(APP_DIR)/junctionFinder
sequenceQuality: $(APP_DIR)/sequenceQuality
wcdInterrogate: $(APP_DIR)/wcdInterrogate

$(APP_DIR)/telomereFinder: build-dirs $(TF_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(TF_OBJ) $(LDFLAGS)

$(APP_DIR)/junctionFinder: build-dirs $(JF_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(JF_OBJ) $(LDFLAGS)

$(APP_DIR)/sequenceQuality: build-dirs $(SQ_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(SQ_OBJ) $(LDFLAGS)

$(APP_DIR)/wcdInterrogate: build-dirs $(WCD_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(WCD_OBJ) $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

.PHONY: all build build-dirs clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

build-dirs:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(OBJ_DIR)/core
	@mkdir -p $(OBJ_DIR)/telomereFinder
	@mkdir -p $(OBJ_DIR)/junctionFinder

debug: CXXFLAGS += -DDEBUG -g
debug: clean all

release: CXXFLAGS += -O2
release: clean all

clean:
	-@rm -rf $(BUILD)
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
