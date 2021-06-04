#include <cstdlib>
#include <cmath>
#include <ctime>
#include <QVector>

QVector<double> n;

double generateWhiteNoiseSample();

QVector<double> generateWhiteNoise(int size) {
    srand(time(nullptr));
    n.resize(size);
    for (int index = 0; index < size; index++)
        n[index] = generateWhiteNoiseSample();

    return n;
}

double generateWhiteNoiseSample() {
    double tempValue;
    double average;
    double result;
    bool p = true;

    while (p) {
        average = (rand() / ((double) RAND_MAX));

        if (average == 0) {
            p = true;
        } else {
            p = false;
        }
    }

    tempValue = cos((2.0 * (double) M_PI) * rand() / ((double) RAND_MAX));
    result = sqrt(-2.0 * log(average)) * tempValue;

    return result;
}
