TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    crystal.cpp \
    atom.cpp \
    verletalgo.cpp \
    zigrandom.c \
    zignor.c

HEADERS += \
    crystal.h \
    atom.h \
    verletalgo.h \
    zigrandom.h \
    zignor.h

LIBS += -larmadillo
