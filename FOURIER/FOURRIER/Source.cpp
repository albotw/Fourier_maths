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
void printVector(vector<complex<double>>* data)
{
    vector<complex<double>>::iterator it;
    for (it = data->begin(); it != data->end(); ++it)
    {
        cout << *it << endl;
    }
}

/*fait un effet miroir pour les bits: 100 -> 001, 101 -> 101*/
/*chiffre à inverser, nombre de bits du chiffre*/
unsigned long mirror(unsigned long data, unsigned long nBits)
{
    unsigned long output = 0l;
    for (int i = 0; i < nBits; i++)
    {
        if (data & (1 << i))
            output |= 1 << ((nBits - 1) - i);
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
complex<double> DFT(complex<double> c1, complex<double> c2, int k, int n, bool isEven, bool invert)
{
    complex<double> omega = 0;
    if (invert)
        omega = exp((-2 * M_PI * I * (double)k) / (double)n);
    else 
        omega = exp((2 * M_PI * I * (double)k) / (double)n);

    if (isEven)
        return c1 + (omega * c2);
    else
        return c1 - (omega * c2);
}

/*
    On calcule le nombre de cycles total à partir du nombre d'échantillons

    pour chaque cycle (1 à cyclesTotal inclus)
        on copie out dans in
        on met la largeur papillon comme étant cycle²
        pour chaque DTF à faire (de 0 à nSamples)
            Si le compteur de DTF est < la moitié de la largeur du papillon:
                out[i] = DTF(in[i], in[i + offset], compteurPap, largeurPap)
            Si le compteur de DTF est > que la moitié de la largeur du papillon:
                out[i] = DTF(in[i - offset], in[i], compteurPap, largeurPap)
            si compteurPap++ >= largeurPap
                compteurPap = 0
*/

vector<complex<double>>* cooleyTurkey(vector<complex<double>> * data, bool inverted)
{

    unsigned int samples = data->size();
    unsigned int nCycles = getNumberCycles(samples);
    vector<complex<double>>* in = new vector<complex<double>>();
    vector<complex<double>>* out = new vector<complex<double>>();
    unsigned int maxExponent = 1;
    unsigned int offset = 1;
    
    for (int i = 0; i < data->size(); i++)
    {
        in->push_back(data->at(mirror(i, 3)));
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

        for (int ftd = 0; ftd < samples; ftd++)
        {
            if (butterflyCounter < (butterflyWidth / 2))
                out->push_back(DFT(in->at(ftd), in->at(ftd + offset), exponent, butterflyWidth, true, inverted));
            else if (butterflyCounter >= (butterflyWidth / 2))
                out->push_back(DFT(in->at(ftd - offset), in->at(ftd), exponent, butterflyWidth, false, inverted));

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

int main()
{
    vector<complex<double>> test1(8);
    test1[0] = complex<double>(0.0, 0.0);
    test1[1] = complex<double>(1.0, 0.0);
    test1[2] = complex<double>(0.0, 0.0);
    test1[3] = complex<double>(0.0, 0.0);
    test1[4] = complex<double>(0.0, 0.0);
    test1[5] = complex<double>(0.0, 0.0);
    test1[6] = complex<double>(0.0, 0.0);
    test1[7] = complex<double>(0.0, 0.0);

    vector<complex<double>>* result1 = cooleyTurkey(&test1, false);
    printVector(result1);

    cout << endl;

    vector<complex<double>>* result2 = cooleyTurkey(result1, true);
    printVector(result2);
}