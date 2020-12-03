#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <complex>
#include <vector>

using namespace std;

const complex<double> I(0.0, 1.0);

/*Afficher un chiffre en binaire*/
/*chiffre à afficher, nombre de bits, caractère pour 0, caractère pour 1*/
void affBin(unsigned long data, unsigned char nBits, char car0, char car1)
{
    for (int i = nBits - 1; i >= 0; i--)
    {
        if (!(data & (1 << i)))
        {
            cout << car0;
        }
        else
        {
            cout << car1;
        }
    }
    cout << endl;
}

/*Affiche un vecteur de complexes*/
void printVector(vector<complex<double>>* data, int M, int N)
{
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << data->at(i + M * j) << " ";
        }
        cout << endl;
    }
}

/*fait un effet miroir pour les bits: 100 -> 001, 101 -> 101*/
/*chiffre à inverser, valeur maximale du chiffre*/
unsigned int mirror(unsigned int data, unsigned int max)
{
    /*comme on sait que max est samples, c'est forcément une puissance de deux*/
    /*Aussi, elle représente le nombre de bits pour la valeur maximale, évitant de modifier les 16 bits et de faire des dépassements de buffer.*/
    unsigned int bits = 0;
    while ((max & (1 << bits)) == 0)
    {
        bits++;
    }

    unsigned int output = 0;
    for (int i = 0; i < bits; i++)
    {
        if (data & (1 << i))
            output |= 1 << ((bits - 1) - i);
    }

    return output;
}

/*Pour compter combien de cycles de recomposition on doit effectuer (2^X = nSamples)*/
unsigned int getNumberCycles(unsigned int nSamples)
{
    unsigned int i = 0;
    while ((nSamples & (1 << i)) == 0)
        i++;

    return i;
}

/*Mini DFT entre deux valeurs*/
/*1er chiffre complexe, 2e chiffre complexe, exposant, indice, traitement pair|impair, DFT normale|inversée*/
complex<double> MiniTFD(complex<double> c1, complex<double> c2, int k, int n, bool isEven, bool invert)
{
    complex<double> omega = 0;
    if (invert)
        omega = exp((2 * M_PI * I * (double)k) / (double)n);
    else 
        omega = exp((-2 * M_PI * I * (double)k) / (double)n);

    if (isEven)
        return c1 + (omega * c2);
    else
        return c1 - (omega * c2);
}

/*
* On met l'exposant max à 1
* On met le décalage du papillon par rapport à in à 1.
* On récupère le nombre d'échantillon
* On calcule le nombre de cycles.
* On copie data dans in avec l'inversion des bits
* 
* Pour chaque cycle (de 1 à nCycle inclus)
*   Si on est pas sur le premier cycle
*       on vide in
*       on copie out dans in
*       on vide out.
* 
*   On met la largeur du papillon comme étant 2^cycle
*   On met le compteur du butterfly à 0.
*   On met le compteur de l'exposant à 0.
*/

vector<complex<double>>* FFT_1D(vector<complex<double>> * data, bool inverted)
{
    unsigned int samples = data->size();
    unsigned int nCycles = getNumberCycles(samples);
    vector<complex<double>>* in = new vector<complex<double>>();
    vector<complex<double>>* out = new vector<complex<double>>();
    unsigned int maxExponent = 1;
    unsigned int offset = 1;
    
    for (int i = 0; i < samples; i++)
    {
        in->push_back(data->at(mirror(i, samples)));
    }
    
    for (int cycle = 1; cycle <= nCycles; cycle++)
    {
        if (cycle != 1)
        {
            in->clear();
            *in = *out;
            out->clear();
        }

        unsigned int butterflyWidth = pow(2.0, cycle);
        unsigned int butterflyCounter = 0;
        unsigned int exponent = 0;

        for (int tfd = 0; tfd < samples; tfd++)
        {
            if (butterflyCounter < (butterflyWidth / 2))
                out->push_back(MiniTFD(in->at(tfd), in->at(tfd + offset), exponent, butterflyWidth, true, inverted));
            else if (butterflyCounter >= (butterflyWidth / 2))
                out->push_back(MiniTFD(in->at(tfd - offset), in->at(tfd), exponent, butterflyWidth, false, inverted));

            if (++butterflyCounter >= butterflyWidth)
                butterflyCounter = 0;

            if (++exponent >= maxExponent)
                exponent = 0;
            
        }
        maxExponent *= 2;
        offset *= 2;
    }

    if (inverted)
    {
        for (int i = 0; i < out->size(); i++)
        {
            out->at(i) = out->at(i) / double(samples);
        }
    }
    
    delete in;
    return out;
}

vector<complex<double>>* FFT_2D(vector<complex<double>>* data, int M, int N, bool invert)
{
    vector<complex<double>>* result = new vector<complex<double>>();
    
    /*Première série de FFT_1D sur les lignes*/
    for (int i = 0; i < N; i++)
    {
        vector<complex<double>>* temp = new vector<complex<double>>();
        temp->insert(temp->begin(), data->begin() + i * M, data->begin() + i * M + M);
        temp = FFT_1D(temp, invert);
        result->insert(result->end(), temp->begin(), temp->end());
        delete temp;
    }

    vector<complex<double>>* output = new vector<complex<double>>(M * N);
    /*Seconde série de FFT_1D sur les colonnes*/
    for (int i = 0; i < M; i++)
    {
        vector<complex<double>>* temp = new vector<complex<double>>();
        //on prend la colonne i de chaque ligne j
        for (int j = 0; j < N; j++)
        {
            temp->push_back(result->at(i + j * M));
        }

        temp = FFT_1D(temp, invert);

        for (int j = 0; j < N; j++)
        {
            output->at(i + j * M) = temp->at(j);
        }
        delete temp;

    }
    delete result;
    return output;
}

vector<complex<double>> * TFD_2D(vector<complex<double>> * data, int M, int N, bool invert)
{
    vector<complex<double>>* out = new vector<complex<double>>();

    for (int u = 0; u < M; u++)
    {
        for (int v = 0; v < M; v++)
        {
            complex<double> result = (0, 0);
            for (int x = 0; x < M; x++)
            {
                for (int y = 0; y < N; y++)
                {
                    double n1 = (double)x * (double)u;
                    double n2 = (double)y * (double)v;
                    double coeff = n1 / M + n2 / N;
                    if (!invert)
                        result += data->at(y + M * x) * exp((-2 * M_PI * I * coeff));
                    else
                        result += data->at(y + M * x) * exp((2 * M_PI * I * coeff));
                }
            }
            out->push_back(result);
        }
    }

    if (invert)
    {
        for (int i = 0; i < M * N; i++)
        {
            out->at(i) = out->at(i) / ((double)M * (double)N);
        }
    }

    return out;
}

vector<complex<double>>* TFD_1D(vector<complex<double>>* data)
{
    int size = data->size();
    vector<complex<double>>* out = new vector<complex<double>>();
    for (int i_out = 0; i_out < size; i_out++)
    {
        complex<double> result = (0,0);
        for (int i = 0; i < size; i++)
        {
            result += data->at(i) * exp((-2 * M_PI * I * ((double)i *(double)i_out)) / (double)size);
        }
        out->push_back(result);
    }

    return out;
}

int main()
{
    vector<complex<double>> test1(4);
    test1[0] = complex<double>(0.0, 0.0);
    test1[1] = complex<double>(1.0, 0.0);
    test1[2] = complex<double>(0.0, 0.0);
    test1[3] = complex<double>(0.0, 0.0);
    cout << "Signal d'origine:" << endl;
    printVector(&test1, 1, 4);
    cout << endl;

    vector<complex<double>>* result1 = FFT_1D(&test1, false);
    cout << "FFT_1D: " << endl;
    printVector(result1, 1, 4);
    cout << endl;

    vector<complex<double>>* result2 = FFT_1D(result1, true);
    cout << "FFT_1D inverse: " << endl;
    printVector(result2, 1, 4);
    cout << endl;

    vector<complex<double>> test2(4);
    test2[0] = complex<double>(1.0, 0.0);
    test2[1] = complex<double>(1.0, 0.0);
    test2[2] = complex<double>(0.0, 0.0);
    test2[3] = complex<double>(0.0, 0.0);
    cout << "Matrice d'origine:" << endl;
    printVector(&test2, 2, 2);
    cout << endl;

    vector<complex<double>>* result3 = TFD_2D(&test2, 2, 2, false);
    cout << "TFD_2D: " << endl;
    printVector(result3, 2, 2);
    cout << endl;

    vector<complex<double>>* result4 = FFT_2D(&test2, 2, 2, false);
    cout << "FFT_2D: " << endl;
    printVector(result4, 2, 2);
    cout << endl;

    vector<complex<double>>* result5 = FFT_2D(result4, 2, 2, true);
    cout << "FFT_2D inverse: " << endl;
    printVector(result5, 2, 2);
    cout << endl;

}