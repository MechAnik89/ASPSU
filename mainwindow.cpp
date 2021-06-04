#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "white_noise.h"
#include "DFT.h"
#include <QtCore>
#include <QtGui>
#include <QPushButton>
#include <QLineEdit>
#include "main.cpp"
#include <iostream>
#include <time.h>
#include <math.h>
#include <QVector>
#include <QPen>
#include <random>

using namespace std;

void generate_signals();
QVector<double> work_of_system(QVector<double> sig);
void show_signals();
double SNR(QVector<double> signal, QVector<double> error);
void myDFT(QVector<double> src);

QVector<double> signal1(1), signal2(1), signal3(1), summary_signal(1), t(1), dft(1), dft_sin(1), f_Four(1), signal8(1), signal16(1);
QVector<double> noise(1), err8(1), err16(1), err8_aft(1), err16_aft(1);
QVector<double> Y_signal, Y_dpf, Y8(1), Y16(1), Y_noise(1);
QVector<double> Y_divided_to_X;
QVector<double> dft_mod(1), dft_arg(1), freqs(1);
double PSD;

int Fs=4000;
float max_time=0.01;
float A1=100, A2=100, A3=100;
int f1=800, f2=1000, f3=2000;
int count_of_samples=(int)(max_time*Fs);

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    srand(time(NULL));
    ui->FreqDiscrLineEdit->setText(QString::number(Fs));
    ui->F1lineEdit->setText(QString::number(f1));
    ui->F2lineEdit->setText(QString::number(f2));
    ui->F3lineEdit->setText(QString::number(f3));
    ui->MaxTimelineEdit->setText(QString::number(max_time));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_GenerateButton_clicked()
{
    Fs=ui->FreqDiscrLineEdit->text().toInt();
    f1=ui->F1lineEdit->text().toInt();
    f2=ui->F2lineEdit->text().toInt();
    f3=ui->F3lineEdit->text().toInt();
    max_time=ui->MaxTimelineEdit->text().toFloat();
    count_of_samples=(int)(max_time*Fs);

    dft.resize(count_of_samples);
    dft_sin.resize(count_of_samples);
    Y_divided_to_X.resize(count_of_samples);
    Y_signal.resize(count_of_samples);
    Y_dpf.resize(count_of_samples);
    dft_arg.resize(count_of_samples);
    dft_mod.resize(count_of_samples);
    freqs.resize(count_of_samples/2);
    Y_noise.resize(count_of_samples);

    std::cout<<Fs<<std::endl;
    std::cout<<max_time<<std::endl;
    std::cout<<count_of_samples<<std::endl;

    std::cout<<"Fs changed to "<<Fs<<std::endl<<"f1 changed to "<<f1<<std::endl<<"f2 changed to "<<f2<<std::endl<<"f3 changed to "<<f3<<std::endl<<std::endl;

    generate_signals();
    std::cout<<"First signal: lenght=" << signal1.length() <<std::endl;
    std::cout<<"Second signal: lenght=" << signal2.length() <<std::endl;
    std::cout<<"Third signal: lenght=" << signal3.length() <<std::endl;
    std::cout<<"Noise: lenght=" << noise.length() <<std::endl;

    //dft_sin=DFT(summary_signal);
    //dft=DFT(noise);
    f_Four.resize(dft.length());
    for(int i=0; i<dft.length(); i++)
        f_Four[i]=(double)i*100;
    Y_signal=work_of_system(summary_signal);
    Y8=work_of_system(signal8);
    Y16=work_of_system(signal16);
    Y_noise=work_of_system(noise);

    //Y_dpf=DFT(Y_signal);
    //myDFT(summary_signal);
    //myDFT(Y_signal);
    //myDFT(signal8);
    //myDFT(Y16);
    myDFT(noise);

    for(int i=0; i<(count_of_samples/2); i++){
        freqs[i]=i*Fs/count_of_samples;
    }

    double snr8=SNR(summary_signal, err8);
    double snr16=SNR(summary_signal, err16);

    err8_aft.resize(count_of_samples);
    err16_aft.resize(count_of_samples);

    for(int i=0; i<count_of_samples; i++){
        err8_aft[i]=Y_signal[i]-Y8[i];
        err16_aft[i]=Y_signal[i]-Y16[i];
    }

    double snr8_aft=SNR(Y_signal, err8_aft);
    double snr16_aft=SNR(Y_signal, err16_aft);

    cout<<endl<<"snr8="<<snr8<<endl<<"snr16="<<snr16<<endl;
    cout<<endl<<"snr8_aft="<<snr8_aft<<endl<<"snr16_aft="<<snr16_aft<<endl<<endl;

    ui->widgetGraph1->clearGraphs();
    ui->widgetGraph1->addGraph();
    ui->widgetGraph1->graph(0)->setPen(QPen(Qt::blue));
    ui->widgetGraph1->graph(0)->setData(t, noise);
    //ui->widgetGraph1->graph(0)->setData(t, signal1);
    //ui->widgetGraph1->graph(0)->setData(t, summary_signal);
    //ui->widgetGraph1->graph(0)->setData(t, signal16);
    //ui->widgetGraph1->graph(0)->setData(f_Four, noise);
    //ui->widgetGraph1->graph(0)->setData(freqs, dft_mod);
    //ui->widgetGraph1->graph(0)->setData(t, Y_signal);
    //ui->widgetGraph1->graph(0)->setData(t, err16);
    ui->widgetGraph1->xAxis->setLabel("t");
    ui->widgetGraph1->yAxis->setLabel("y");
    //ui->widgetGraph1->xAxis->setRange(0, max_time);
    //ui->widgetGraph1->yAxis->setRange(-A1*1, A1*1);
    ui->widgetGraph1->yAxis->setAutoTickStep(true);
    ui->widgetGraph1->graph(0)->rescaleAxes();
    ui->widgetGraph1->replot();

    ui->widgetGraph2->clearGraphs();
    ui->widgetGraph2->addGraph();
    ui->widgetGraph2->graph(0)->setPen(QPen(Qt::blue));
    ui->widgetGraph2->graph(0)->setData(freqs, dft_mod);
    //ui->widgetGraph2->graph(0)->setData(freqs, dft_arg);
    //ui->widgetGraph2->graph(0)->setData(t, signal2);
    //ui->widgetGraph2->graph(0)->setData(t, err16_aft);
    //ui->widgetGraph2->graph(0)->setData(t, Y_signal);
    //ui->widgetGraph2->graph(0)->setData(freqs, dft_mod);
    //ui->widgetGraph2->graph(0)->setData(t, Y16);
    ui->widgetGraph2->xAxis->setLabel("f");
    ui->widgetGraph2->yAxis->setLabel("y");
    //ui->widgetGraph2->xAxis->setRange(0, max_time);
    //ui->widgetGraph2->yAxis->setRange(-A2*1, A2*1);
    ui->widgetGraph2->yAxis->setAutoTickStep(true);
    ui->widgetGraph2->graph(0)->rescaleAxes();
    ui->widgetGraph2->replot();

    ui->widgetGraph3->clearGraphs();
    ui->widgetGraph3->addGraph();
    ui->widgetGraph3->graph(0)->setPen(QPen(Qt::blue));
    //ui->widgetGraph3->graph(0)->setData(t, signal3);
    ui->widgetGraph3->graph(0)->setData(t, Y_noise);
    //ui->widgetGraph3->graph(0)->setData(f_Four, err16);
    ui->widgetGraph3->xAxis->setLabel("t");
    ui->widgetGraph3->yAxis->setLabel("y");
    ui->widgetGraph3->xAxis->setRange(0, max_time);
    //ui->widgetGraph3->yAxis->setRange(-A3*1, A3*1);
    ui->widgetGraph3->yAxis->setAutoTickStep(true);
    ui->widgetGraph3->graph(0)->rescaleAxes();
    ui->widgetGraph3->replot();

    myDFT(Y_noise);

    ui->widgetGraphSum->clearGraphs();
    ui->widgetGraphSum->addGraph();
    ui->widgetGraphSum->graph(0)->setPen(QPen(Qt::blue));
    //ui->widgetGraphSum->graph(0)->setData(t, summary_signal);
    ui->widgetGraphSum->graph(0)->setData(freqs, dft_mod);
    ui->widgetGraphSum->xAxis->setLabel("f");
    ui->widgetGraphSum->yAxis->setLabel("y");
    ui->widgetGraphSum->xAxis->setRange(0, max_time);
    //ui->widgetGraphSum->yAxis->setRange(-(A1+A2+A3)*1, (A1+A2+A3)*1);
    ui->widgetGraphSum->yAxis->setAutoTickStep(true);
    ui->widgetGraphSum->graph(0)->rescaleAxes();
    ui->widgetGraphSum->replot();

    ui->widgetGraphDFT->clearGraphs();
    ui->widgetGraphDFT->addGraph();
    ui->widgetGraphDFT->graph(0)->setPen(QPen(Qt::blue));
    //ui->widgetGraphDFT->graph(0)->setData(t, Y_signal);
    ui->widgetGraphDFT->xAxis->setLabel("t");
    ui->widgetGraphDFT->yAxis->setLabel("y");
    ui->widgetGraphDFT->xAxis->setRange(0, max_time);
    //ui->widgetGraphDFT->yAxis->setRange(-A1*1, A1*1);
    ui->widgetGraphDFT->yAxis->setAutoTickStep(true);
    ui->widgetGraphDFT->graph(0)->rescaleAxes();
    ui->widgetGraphDFT->replot();
}

void generate_signals(){
    double max=0;
    t.resize(count_of_samples);
    signal1.resize(count_of_samples);
    signal2.resize(count_of_samples);
    signal3.resize(count_of_samples);

    signal8.resize(count_of_samples);
    signal16.resize(count_of_samples);
    err8.resize(count_of_samples);
    err16.resize(count_of_samples);

    summary_signal.resize(count_of_samples);
    noise.resize(count_of_samples);
    noise=generateWhiteNoise(noise.length());

    for(int i=0; i<count_of_samples; i++){
        t[i]=i*(1.0/Fs);
        signal1[i]=sin(2*M_PI*f1*t[i]);
        signal2[i]=sin(2*M_PI*f2*t[i]);
        signal3[i]=sin(2*M_PI*f3*t[i]);
        summary_signal[i]=(signal1[i]+signal2[i]+signal3[i])/3;
        signal8[i]=trunc(summary_signal[i]*(pow(2, 8)-1))/(pow(2, 8)-1);
        signal16[i]=trunc(summary_signal[i]*(pow(2, 16)-1))/(pow(2, 16)-1);
        err8[i]=summary_signal[i]-signal8[i];
        err16[i]=summary_signal[i]-signal16[i];

        if(fabs(noise[i])>max)
            max=fabs(noise[i]);
    }
    for(int i=0; i<count_of_samples; i++){
        noise[i]/=max;
    }
}

void show_signals(){
    std::cout<<"First signal: lenght=" << signal1.length() <<std::endl;
    for(int i=0; i<count_of_samples; i++){
        std::cout<<signal1[i]<<std::endl;
    }
    std::cout<<std::endl;
}

QVector<double> work_of_system(QVector<double> sig){
    QVector<double> res;
    res.resize(count_of_samples);
    double y1=0, y2=0, y1_1=0, y1_2=0, y2_1=0, y2_2=0, x_1=0, x_2=0;
    for(int i=0; i<count_of_samples; i++){
        y2_2=y2_1;
        y1_2=y1_1;
        y2_1=y2;
        y1_1=y1;

        y1=0.6293*sig[i]-1.2586*x_1+0.6293*x_2+0.9747*y1_1-0.5426*y1_2;
        y2=0.4754*y1-0.9508*y1_1+0.4754*y1_2+0.7363*y2_1-0.1654*y2_2;
        res[i]=y2;
    }
    return res;
}

double SNR(QVector<double> signal, QVector<double> error){
    double snr=0, sum1=0, sum2=0;
    for(int i=0; i<signal.size(); i++){
        sum1+=signal[i]*signal[i];
        sum2+=error[i]*error[i];
    }
    sum1/=signal.size();
    sum2/=signal.size();
    snr=10*log(sum1/sum2);

    return snr;
}

void myDFT(QVector<double> src){
    double abs, arg;
    for(int k=0; k<count_of_samples; k++){
        abs=0;
        arg=0;
        for(int n=0; n<count_of_samples; n++){
            arg+=cos(2*M_PI*(double)n*k/count_of_samples)*src[n];
            abs+=sin(2*M_PI*(double)n*k/count_of_samples)*src[n]*(-1.0);
        }
        dft_arg[k]=arg;
        dft_mod[k]=fabs(abs);
    }
    for(int i=0; i<(count_of_samples/2); i++){
        PSD+=dft_mod[i]*dft_mod[i];
    }
    PSD/=count_of_samples/2;
    std::cout<<"PSD="<<PSD<<std::endl;
}
