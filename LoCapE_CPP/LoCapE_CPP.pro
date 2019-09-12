QT += core
QT -= gui

CONFIG += c++11

TARGET = testEigenProcess
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    locape.cpp \

# The following define makes your compiler emit warnings if you use
# any feature of Qt which as been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Eigen
INCLUDEPATH += \
/media/jianboma/Jianbo_Research/NMR_project/LoCapE_CPP/LoCapE_CPP
# openMP
QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp
# intel MKL
unix:INCLUDEPATH += /opt/intel/mkl/include
unix:LIBS += -L/opt/intel/mkl/lib/intel64 \
    -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
    -L/home/jianboma/intel/mkl/lib/intel64 \
    -liomp5 -lpthread -ldl -lm

HEADERS += \
    locape.h
