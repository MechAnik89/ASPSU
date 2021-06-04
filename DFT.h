#include <math.h>
#include <QVector>

QVector<double> mas_dft;

QVector<double> DFT(QVector<double> arr) {
    mas_dft.resize(arr.length());
    double Re = 0, Im = 0, Arg = 0;
    for (int i = 0; i < arr.length(); i++) {
    Re = 0;
    Im = 0;
    for (int j = 0; j < arr.length(); j++) {
    Arg = 2 * M_PI * i * j / arr.length();
    Re += cos(Arg) * arr[j];
    Im -= sin(Arg) * arr[j];
    }
    mas_dft[i] = (sqrt(pow(Re, 2) + pow(Im, 2))) * 2 / arr.length();
    }
    return mas_dft;
}

