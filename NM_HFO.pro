TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

TARGET = HFO.cpp

SOURCES +=  CA3_Column.cpp	\
	    Cortical_Column.cpp \
	    HFO_mex.cpp		\
	    HFO.cpp	    

HEADERS +=  CA3_Column.h	\
	    Cortical_Column.h	\
	    Data_Storage.h	\
	    ODE.h		\
	    Random_Stream.h
    

QMAKE_CXXFLAGS += -std=c++11 -O3

SOURCES -= HFO_mex.cpp
