TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \

HEADERS += \
    common.h \
    operators.h \
    tests.h \
    tmatrix.h \
    tquaternion.h \
    tsymmetricmatrix.h \
    tvector.h
