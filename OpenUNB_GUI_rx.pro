#-------------------------------------------------
#
# Project created by QtCreator 2020-11-10T10:25:37
#
#-------------------------------------------------

QT       += core gui charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = OpenUNB_GUI
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++17

SOURCES += \
        main.cpp \
        mainwindow.cpp \
        SDRDevInfo.cpp \
        SoapyEnum.cpp \
        demodulator.cpp \
        singenerator.cpp \
        pcscl_noperm.cpp \
        crc_ok.cpp \
    udp_reciever.cpp \
    preamblepoint.cpp \
    decoder.cpp \
    VectorHelper.cpp

HEADERS += \
        mainwindow.h \
        SDRDevInfo.h \
        SoapyEnum.h \
        demodulator.h \
        fdacoefs_125K_to_100.h \
        fdacoefs_1M_to_125K.h \
        fdacoefs_1M_to_100.h \
        singenerator.h \
    pcscl_noperm.h \
    crc_ok.h \
    udp_reciever.h \
    preamblepoint.h \
    decoder.h \
    VectorHelper.h

FORMS += \
        mainwindow.ui

INCLUDEPATH +=  /usr/local/include

LIBS += -L/usr/local/lib \
        -L/usr/lib \
        -lfftw3f \
        -lSoapySDR \
        -pthread

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
