# ----------------------------------------------------
# This file is generated by the Qt Visual Studio Add-in.
# ------------------------------------------------------

TEMPLATE = lib
TARGET = UEBHydroCoupleComponent
DESTDIR = ../Mac_x86_64/Debug
QT += core
CONFIG += debug
DEFINES += WIN64 QT_DLL UEBHYDROCOUPLECOMPONENT_LIB
INCLUDEPATH += ./GeneratedFiles \
    . \
    ./GeneratedFiles/Debug
PRECOMPILED_HEADER = stdafx.h
DEPENDPATH += .
MOC_DIR += ./GeneratedFiles/debug
OBJECTS_DIR += debug
UI_DIR += ./GeneratedFiles
RCC_DIR += ./GeneratedFiles
include(UEBHydroCoupleComponent.pri)



win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../../../../usr/local/Cellar/netcdf/4.3.3.1/lib/release/ -lnetcdf
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../../../../usr/local/Cellar/netcdf/4.3.3.1/lib/debug/ -lnetcdf
else:unix: LIBS += -L$$PWD/../../../../../../../../usr/local/Cellar/netcdf/4.3.3.1/lib/ -lnetcdf

INCLUDEPATH += $$PWD/../../../../../../../../usr/local/Cellar/netcdf/4.3.3.1/include
DEPENDPATH += $$PWD/../../../../../../../../usr/local/Cellar/netcdf/4.3.3.1/include



win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/lib/release/ -lmpicxx
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/lib/debug/ -lmpicxx
else:unix: LIBS += -L$$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/lib/ -lmpicxx

INCLUDEPATH += $$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/include
DEPENDPATH += $$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/lib/release/ -lmpi
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/lib/debug/ -lmpi
else:unix: LIBS += -L$$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/lib/ -lmpi

INCLUDEPATH += $$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/include
DEPENDPATH += $$PWD/../../../../../../../../usr/local/Cellar/mpich2/3.1.4_1/include