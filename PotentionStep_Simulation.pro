TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS += -fopenmp
LIBS += -lgomp -lpthread

SOURCES += main.cpp \
    functions.cpp \
    coefficient_thomas.cpp \
    iodata.cpp \
    simdata.cpp

HEADERS += \
    coefficient_thomas.h \
    iodata.h \
    simdata.h \
    functions.h
