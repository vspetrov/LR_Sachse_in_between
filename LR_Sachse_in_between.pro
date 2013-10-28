TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    LR_lattice.cpp \
    LR_cell.cpp \
    Sachse_fibrolast.cpp

HEADERS += \
    LR_lattice.h \
    LR_cell.h \
    Sachse_fibroblast.h

QMAKE_CXXFLAGS = -std=c++0x
QMAKE_CC = mpicc
QMAKE_CXX = mpic++
QMAKE_LINK = mpic++
