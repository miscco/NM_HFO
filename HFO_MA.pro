TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES +=  Main.cpp \
            HFO.cpp \
            H_Column.cpp \
            C_Column.cpp

HEADERS +=  \
            ODE.h \
            saves.h \
            Stimulation.h \
            H_Column.h \
            C_Column.h

QMAKE_CXXFLAGS += -std=c++11 -O3

SOURCES -= HFO.cpp
