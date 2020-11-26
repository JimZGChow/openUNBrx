#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow) {
    ui->setupUi(this);

    connect(ui->pushButton_connect, SIGNAL (released()), this, SLOT (evtConnetc()));
    connect(ui->checkBox_run_1mhz, SIGNAL(stateChanged(int)), this, SLOT(evtCheckBoxRun(int)));
    connect(ui->checkBox_run_100hz, SIGNAL(stateChanged(int)), this, SLOT(evtCheckBoxRun(int)));

    connect(ui->checkBox_100_freq, SIGNAL(stateChanged(int)), this, SLOT(evtCheckBoxIQFP(int)));
    connect(ui->checkBox_100_i, SIGNAL(stateChanged(int)), this, SLOT(evtCheckBoxIQFP(int)));
    connect(ui->checkBox_100_q, SIGNAL(stateChanged(int)), this, SLOT(evtCheckBoxIQFP(int)));
    connect(ui->checkBox_100_phase, SIGNAL(stateChanged(int)), this, SLOT(evtCheckBoxIQFP(int)));

    connect(ui->lineEdit_100hz, SIGNAL(returnPressed()), this, SLOT(evtChangeChannel()));
    connect(ui->lineEdit_freq, SIGNAL(returnPressed()), this, SLOT(evtChangeFreq()));
    connect(ui->lineEdit_div, SIGNAL(returnPressed()), this, SLOT(evtChangeDiv()));
    connect(ui->pushButton_save, SIGNAL(released()), this, SLOT(evtSaveToFile()));

    qChart1M = new QChart();
    qLineSeries1M = new QLineSeries();
    qChart100 = new QChart();
    qLineSeries100_I = new QLineSeries();
    qLineSeries100_Q = new QLineSeries();
    qLineSeries100_Phase = new QLineSeries();
    qLineSeries100_Freq = new QLineSeries();
    qLineSeries100_En = new QLineSeries();
    qLineSeries100_Corr = new QLineSeries();
    qLineSeries100_Bit = new QLineSeries();
    recvDataArray = new QByteArray();
    QChartView *chartView = new QChartView(qChart1M);
    QChartView *chartView2 = new QChartView(qChart100);

    qChart1M->addSeries(qLineSeries1M);
    QValueAxis *axisX = new QValueAxis;
    axisX->setRange(-CHAN_NUM/2, CHAN_NUM/2);
    axisX->setLabelFormat("%g");
    axisX->setTitleText("Freq");
    QValueAxis *axisY = new QValueAxis;
    axisY->setRange(-100.00, 20.00);
    axisY->setTitleText("Audio level");
    qChart1M->addAxis(axisX, Qt::AlignBottom);
    qLineSeries1M->attachAxis(axisX);
    qChart1M->addAxis(axisY, Qt::AlignLeft);
    qLineSeries1M->attachAxis(axisY);
    qChart1M->legend()->hide();
    qChart1M->setTitle("Data");

    ui->verticalLayout_1mhz->addWidget(chartView);

    qChart100->addSeries(qLineSeries100_I);
    qChart100->addSeries(qLineSeries100_Q);
    qChart100->addSeries(qLineSeries100_Freq);
    qChart100->addSeries(qLineSeries100_Phase);
    qChart100->addSeries(qLineSeries100_En);
    qChart100->addSeries(qLineSeries100_Corr);
    qChart100->addSeries(qLineSeries100_Bit);
    QValueAxis *axisX2 = new QValueAxis;
    axisX2->setRange(0, maxIQ);
    axisX2->setLabelFormat("%g");
    axisX2->setTitleText("Samples");
    axisX2->setTickCount(101);
    //axisX2->setTickCount(11);
    QValueAxis *axisY2 = new QValueAxis;
    axisY2->setRange(-1.00, 1.00);
    axisY2->setTitleText("Audio level");
    qChart100->addAxis(axisX2, Qt::AlignBottom);
    qLineSeries100_I->attachAxis(axisX2);
    qLineSeries100_Q->attachAxis(axisX2);
    qLineSeries100_Phase->attachAxis(axisX2);
    qLineSeries100_Freq->attachAxis(axisX2);
    qLineSeries100_En->attachAxis(axisX2);
    qLineSeries100_Corr->attachAxis(axisX2);
    qLineSeries100_Bit->attachAxis(axisX2);
    qChart100->addAxis(axisY2, Qt::AlignLeft);
    qLineSeries100_I->attachAxis(axisY2);
    qLineSeries100_Q->attachAxis(axisY2);
    qLineSeries100_Phase->attachAxis(axisY2);
    qLineSeries100_Freq->attachAxis(axisY2);
    qLineSeries100_En->attachAxis(axisY2);
    qLineSeries100_Corr->attachAxis(axisY2);
    qLineSeries100_Bit->attachAxis(axisY2);
    QColor cRed, cGreen, cBlue, cColor, cCorr;
    cRed.setRgb(255, 0, 0);
    cGreen.setRgb(0, 255, 0);
    cBlue.setRgb(0, 0, 0);
    cColor.setRgb(120, 120, 120);
    cCorr.setRgb(21,45,2);
    qLineSeries100_Q->setColor(cRed);
    qLineSeries100_Phase->setColor(cGreen);
    qLineSeries100_Freq->setColor(cBlue);
    qLineSeries100_En->setColor(cColor);
    qLineSeries100_Corr->setColor(cCorr);
    qChart100->legend()->hide();
    qChart100->setTitle("Data");

    ui->verticalLayout_100hz->addWidget(chartView2);

    sinG = new SinCosGenerator(100, 0);

    data100hzGUI = new float[maxIQ*10];

#ifdef FFTW3F
    in = (fftwf_complex*) fftwf_malloc (sempNum * sizeof (fftwf_complex));
    out = (fftwf_complex*) fftwf_malloc (sempNum * sizeof (fftwf_complex));
    p = fftwf_plan_dft_1d(sempNum, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
#endif
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::evtConnected() {
    ui->textBrowser_log->append("Connected!");
    ui->pushButton_connect->setText("Connected");
}


void MainWindow::evtCheckBoxRun(int state ) {
    /*
    if (soc->isOpen()) {
        sendData data;
        data.sendDataType = sendData::sendableData;
        data.fields.fieldSendableData.data_100hz = ui->checkBox_run_100hz->checkState() == 2 ? 1 : 0;
        data.fields.fieldSendableData.data_1mhz = ui->checkBox_run_1mhz->checkState() == 2 ? 1 : 0;

        soc->write((char*)&data, sizeof(sendData));

    }
    */
}

void MainWindow::evtCheckBoxIQFP(int state ) {
    if (ui->checkBox_100_i->checkState() == 2)
        qLineSeries100_I->replace(chanDataI);
    else
        qLineSeries100_I->clear();

    if (ui->checkBox_100_q->checkState() == 2)
        qLineSeries100_Q->replace(chanDataQ);
    else
        qLineSeries100_Q->clear();

    if (ui->checkBox_100_phase->checkState() == 2)
        qLineSeries100_Phase->replace(chanDataPhase);
    else
        qLineSeries100_Phase->clear();

    if (ui->checkBox_100_freq->checkState() == 2)
        qLineSeries100_Freq->replace(chanDataFreq);
    else
        qLineSeries100_Freq->clear();

    if (ui->checkBox_100_en->checkState() == 2)
        qLineSeries100_En->replace(chanDataEn);
    else
        qLineSeries100_En->clear();

    if (ui->checkBox_corr->checkState() == 2)
        qLineSeries100_Corr->replace(chanDataCorr);
    else
        qLineSeries100_Corr->clear();

    if (ui->checkBox_bit->checkState() == 2)
        qLineSeries100_Bit->replace(chanDataBit);
    else
        qLineSeries100_Bit->clear();
}

void MainWindow::wheelEvent(QWheelEvent *event) {
    //qDebug() << "evt " << event->angleDelta();
    //qDebug() << ui->tabWidget->currentIndex();

    if (ui->tabWidget->currentIndex() == 2) {
        //qDebug() << "evt " << event->buttons();

        if (event->buttons() == Qt::LeftButton) { /// X
            if (event->angleDelta().y() < 0) {
                maxIQ += 10;
                qChart100->axes()[0]->setRange(0, maxIQ);

            }
            else {
                if (maxIQ - 10 > 10)
                    maxIQ -= 10;
                qChart100->axes()[0]->setRange(0, maxIQ);
            }


            //soc->write((char*)&data, sizeof(sendData));
        } else {

            if (event->angleDelta().y() < 0) { /// Y
                float oldMin = ((QValueAxis*)qChart100->axes()[1])->min();
                float oldMax = ((QValueAxis*)qChart100->axes()[1])->max();
                qChart100->axes()[1]->setRange(oldMin * 1.25, oldMax * 1.25);

            }
            else {
                float oldMin = ((QValueAxis*)qChart100->axes()[1])->min();
                float oldMax = ((QValueAxis*)qChart100->axes()[1])->max();
                qChart100->axes()[1]->setRange(oldMin * 0.75, oldMax * 0.75);
            }
        }
    }
}

void  MainWindow::evtChangeChannel() {
    outChNum = ui->lineEdit_100hz->text().toInt();
    if (outChNum < 0 || outChNum > CHAN_NUM)
        outChNum = 0;


    //soc->write((char*)&data, sizeof(sendData));

    qLineSeries100_I->replace(chanDataI);
    qLineSeries100_Q->replace(chanDataQ);
}

void MainWindow::evtChangeDiv() {
    div = ui->lineEdit_div->text().toInt();
    //sinG->setFreq(fr);
}

void MainWindow::evtChangeFreq() {
    float fr = ui->lineEdit_freq->text().toDouble();
    sinG->setFreq(fr);
    update100HzData();
}

float MainWindow::getSelectedFreq() {
    float fr = ui->lineEdit_freq->text().toDouble();
    sinG->setFreq(fr);
    return  fr;
}

void MainWindow::evtSaveToFile() {
    QFile file("data.complex");
    if (!file.open(QIODevice::WriteOnly)) {
    }
    else {
        file.write((char*)data100hzGUI, maxIQ * sizeof(float) * 2);
    }

    file.close();

}

void MainWindow::push100HzData(float* data_in) {

//    for (int i=0; i<maxIQ; i++) {
//        data100hz.push_back(data_in[i*2 + 0]);
//        data100hz.push_back(data_in[i*2 + 1]);
//    }

    data100hz.push_back(data_in[0]);
    data100hz.push_back(data_in[1]);

    if (ui->checkBox_run_100hz->checkState() == 2) {
        if (data100hz.size() >= maxIQ*2) {
            memcpy(data100hzGUI, data100hz.data(), data100hz.size() * sizeof(float));
            update100HzData();

            data100hz.clear();

        }
    }
    else {
        data100hz.clear();
    }
}

void MainWindow::push1MHzData(fftwf_complex *data_in, int size) {
    if (ui->checkBox_run_1mhz->checkState() != 2)
        return;

    QVector<QPointF> m_buffer;
    QPointF max(0, -10000);


    for (int i=0; i < size; i++) {
        //QPoint tmp(i - CHAN_NUM/2, 10*log10f(data_in[i][0]*data_in[i][0] + data_in[i][1]*data_in[i][1]));
        //m_buffer.push_back(tmp);
        int pos;
        if (i < size/2)
            pos = i + size/2;
        else
            pos = i - size/2;

        m_buffer.push_back(QPointF(i - CHAN_NUM/2, 10*log10f(data_in[pos][0]*data_in[pos][0] + data_in[pos][1]*data_in[pos][1])));
        if (m_buffer[i].y() > max.y() && m_buffer[i].x() != -1 && m_buffer[i].x() != 0 && m_buffer[i].x() != 1) {
            max.setY(m_buffer[i].y());
            max.setX(m_buffer[i].x());
        }
    }

    if ( max.y() > -70)
        ui->textBrowser_log->append("Signal on " + QString::number(max.x()) + " ch");


    qLineSeries1M->replace(m_buffer);
}

int MainWindow::getSelectedChannel() {
    return outChNum;
}

void MainWindow::update100HzData() {
    float* dataF = data100hzGUI;
    chanDataI.clear();
    chanDataQ.clear();
    chanDataPhase.clear();
    chanDataFreq.clear();
    chanDataEn.clear();
    chanDataCorr.clear();
    chanDataBit.clear();

    float phase1, phase2;
    //sinG->reset();
    for (int i=0; i < maxIQ; i++) {
        float sI = 1, sQ = 1;

        if (ui->checkBox_mul->checkState() == 2) {
            sI = sinG->nextCos();
            sQ = sinG->nextSin();
            sinG->nextStep();
        }

        chanDataI.push_back(QPointF(maxIQ - 1 - i, dataF[i*2] * sI));
        chanDataQ.push_back(QPointF(maxIQ - 1 - i, dataF[i*2 + 1] * sQ));
        chanDataEn.push_back(QPointF(maxIQ - 1 - i, 10*log10f(dataF[i*2] * sI * dataF[i*2] * sI + dataF[i*2+1] * sQ * dataF[i*2+1] * sQ)));

            //chanDataCorr.push_back(QPointF(maxIQ - 1 - i, -dataF[(i-1)*2] * dataF[i*2 + 1] + dataF[(i-1)*2 + 1] * dataF[i*2]));

        //if (chanDataEn[i].y() > -60) {
        phase1 = atan2f(dataF[i*2 + 1] * sQ, dataF[i*2] * sI);
        chanDataPhase.push_back(QPointF(maxIQ - 1 - i, phase1/div));
        chanDataFreq.push_back(QPointF(maxIQ - 1 - i, (phase1 - phase2)/div));
        float ff = fmodf(phase1 - phase2, 2 * M_PI);
        if (ff < 0)
            ff += M_PI*2;

        if ((ff > M_PI/2) && (ff < 3*M_PI/2))
            chanDataBit.push_back(QPointF(maxIQ - 1 - i, 0.1));
         else
            chanDataBit.push_back(QPointF(maxIQ - 1 - i, -0.1));
        phase2 = phase1;

        if (i > 0)
            chanDataCorr.push_back(QPointF(maxIQ - 1 - i, dataF[i*2] * sI * dataF[(i-1)*2] * sI + dataF[i*2 + 1] * sQ * dataF[(i-1)*2 + 1] * sQ));
        //}
    }

    if (ui->checkBox_100_i->checkState() == 2)
        qLineSeries100_I->replace(chanDataI);
    else
        qLineSeries100_I->clear();

    if (ui->checkBox_100_q->checkState() == 2)
        qLineSeries100_Q->replace(chanDataQ);
    else
        qLineSeries100_Q->clear();

    if (ui->checkBox_100_phase->checkState() == 2)
        qLineSeries100_Phase->replace(chanDataPhase);
    else
        qLineSeries100_Phase->clear();

    if (ui->checkBox_100_freq->checkState() == 2)
        qLineSeries100_Freq->replace(chanDataFreq);
    else
        qLineSeries100_Freq->clear();

    if (ui->checkBox_100_en->checkState() == 2)
        qLineSeries100_En->replace(chanDataEn);
    else
        qLineSeries100_En->clear();

    if (ui->checkBox_corr->checkState() == 2)
        qLineSeries100_Corr->replace(chanDataCorr);
    else
        qLineSeries100_Corr->clear();

    if (ui->checkBox_bit->checkState() == 2)
        qLineSeries100_Bit->replace(chanDataBit);
    else
        qLineSeries100_Bit->clear();
}
