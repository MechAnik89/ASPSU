#include "mainwindow.h"
#include <QApplication>
#include <string>
#include <QVector>
#include <iostream>

void generate_signals();
float* get_address(std::string str);

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;

    w.show();
    return a.exec();
}
