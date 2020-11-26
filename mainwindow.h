#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QFile>
#include <QtCharts/QXYSeries>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QChart>
#include <QtCharts/QValueAxis>
#include <QWheelEvent>
#include <math.h>
#include "singenerator.h"

#define FFTW3F

#ifdef FFTW3F
#include "fftw3.h"
#endif

#define RECIEVE_ALL_DATA    0
#define CHAN_NUM    (1250*2)

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

QT_CHARTS_USE_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    int getSelectedChannel();
    float getSelectedFreq();
    void push100HzData(float* data_in);
    void push1MHzData(fftwf_complex* data_in, int size);

private:
    Ui::MainWindow* ui;
    QChart* qChart1M;
    QLineSeries* qLineSeries1M;
    QChart* qChart100;
    QLineSeries* qLineSeries100_I;
    QLineSeries* qLineSeries100_Q;
    QLineSeries* qLineSeries100_Phase;
    QLineSeries* qLineSeries100_Freq;
    QLineSeries* qLineSeries100_En;
    QLineSeries* qLineSeries100_Corr;
    QLineSeries* qLineSeries100_Bit;

    QByteArray* recvDataArray;
    bool expData;

    QVector<QPointF> chanDataI;
    QVector<QPointF> chanDataQ;
    QVector<QPointF> chanDataPhase;
    QVector<QPointF> chanDataFreq;
    QVector<QPointF> chanDataEn;
    QVector<QPointF> chanDataCorr;
    QVector<QPointF> chanDataBit;
    QVector<float> data100hz;
    float* data100hzGUI;
    int maxIQ = 100;
    int outChNum = 0;
    int div = 1;

    SinCosGenerator* sinG;
#ifdef FFTW3F
    const int sempNum = 100;
    fftwf_complex* in;
    fftwf_complex* out;
    fftwf_plan p;
#endif

    void wheelEvent(QWheelEvent *event);
    void update100HzData();

private slots:
    //void evtConnetc();
    void evtConnected();
    //void evtReadData();
    void evtSaveToFile();
    void evtCheckBoxRun(int state);
    void evtCheckBoxIQFP(int state);
    void evtChangeChannel();
    void evtChangeDiv();
    void evtChangeFreq();
};
#endif // MAINWINDOW_H
