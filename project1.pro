TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    crystal.cpp \
    atom.cpp \
    verletalgo.cpp \
    printing.cpp \
    zigrandom.c \
    zignor.c \
    cell.cpp \
    verletalgo2.cpp

HEADERS += \
    crystal.h \
    atom.h \
    verletalgo.h \
    printing.h \
    zigrandom.h \
    zignor.h \
    cell.h \
    verletalgo2.h

LIBS += -larmadillo
