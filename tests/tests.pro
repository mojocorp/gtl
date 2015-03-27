TEMPLATE = app
TARGET = tests
CONFIG -= app_bundle

DEPENDPATH += . ../gtl
INCLUDEPATH += . ..

# Input
HEADERS += gtest/gtest.h
SOURCES += gtest/gtest-all.cc

SOURCES += main.cpp \
           testBox2.cpp \
           testBox3.cpp \
           testCircle.cpp \
           testComplex.cpp \
           testMatrix3.cpp \
           testMatrix4.cpp \
           testPlane.cpp \
           testQuat.cpp \
           testRay.cpp \
           testSphere.cpp \
           testVec2.cpp \
           testVec3.cpp \
           testVec4.cpp
