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
            printf("%c", car0);
        }
        else
        {
            printf("%c", car1);
        }
    }
    printf("\n");
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

/*Pour compter combien de cycles de recomposition on doit effectuer*/
unsigned int getNumberCycles(unsigned int nSamples)
{
    unsigned int i = 0;
    while ((nSamples & (1 << i)) == 0)
        i++;

    return i;
}
/*Mini DFT entre deux valeurs*/
/*Tableau de 2 complexes, échantillon actuel, nombre d'échantillons*/
complex<double> DFT(complex<double> c1, complex<double> c2, int k, int n, bool isEven)
{
    complex<double> omega = exp((-2 * M_PI * I * (double)k) / (double)n);
    if (isEven)
        return c1 + (omega * c2);
    else
        return c1 - (omega * c2);
}

complex<double> iDFT(complex<double> c1, complex<double> c2, int k, int n, bool isEven)
{
    complex<double> omega = exp((2 * M_PI * I * (double)k) / (double)n);
    if (isEven)
        return c1 + (omega * c2);
    else
        return c1 - (omega * c2);
}

vector<complex<double>> cooleyTurkey(vector<complex<double>> data)
{
    /*
        On calcule le nombre de cycles total à partir du nombre d'échantillons
        
        pour chaque cycle (1 à cyclesTotal inclus)
            on copie out dans in
            on met la largeur papillon comme étant cycle²
            pour chaque DTF à faire (de 0 à nSamples)
                Si le compteur de DTF est < la moitié de la largeur du papillon:
                    out[i] = DTF(in[i], in[i + cycle], compteurPap, largeurPap)
                Si le compteur de DTF est > que la moitié de la largeur du papillon:
                    out[i] = DTF(in[i - cycle], in[i], compteurPap, largeurPap) 
                si compteurPap++ >= largeurPap
                    compteurPap = 0
    */
    unsigned int samples = data.size();
    vector<complex<double>> in(samples);
    for (int i = 0; i < data.size(); i++)
    {
        in[i] = data[mirror(i, 3)];
    }
    
    vector<complex<double>> out(samples);
    unsigned int maxExponent = 1;
    unsigned int butterflyOffset = 1;
    unsigned int nCycles = getNumberCycles(samples);
    for (int cycle = 1; cycle <= nCycles; cycle++)
    {
        if (cycle != 1)
            in = out;

        unsigned int butterflyWidth = pow(2.0, cycle);
        unsigned int butterflyCounter = 0;
        unsigned int exponent = 0;
        //cout << butterflyWidth << endl;
        for (int ftd = 0; ftd < samples; ftd++)
        {
            cout << "but " << butterflyCounter << endl;
            cout << "expo " << exponent << endl;
            //cout << (ftd > (butterflyWidth / 2)) << endl;
            if (butterflyCounter < (butterflyWidth / 2))
            { 
                cout << "c1 " << ftd << " c2 " << (ftd + cycle) << endl;
                out[ftd] = DFT(in[ftd], in[ftd + butterflyOffset], exponent, butterflyWidth, true);
                cout << "out" << out[ftd] << endl;
            }
            else if (butterflyCounter >= (butterflyWidth / 2))
            { 
                cout << "c1 " << (ftd - cycle) << " c2 " << ftd << endl;
                out[ftd] = DFT(in[ftd - butterflyOffset], in[ftd], exponent, butterflyWidth, false);
                cout << "out:"<< out[ftd] << endl;
            }
            cout << endl;

            butterflyCounter++;
            if (butterflyCounter >= butterflyWidth)
            {
                butterflyCounter = 0;
            }

            exponent++;
            if (exponent >= maxExponent)
            {
                exponent = 0;
            }
            
        }
        cout << "--- CYCLE ---" << endl;

        maxExponent *= 2;

        butterflyOffset *= 2;
    }

    
    
    return out;
}

vector<complex<double>> icooleyTurkey(vector<complex<double>> data)
{
    /*
        On calcule le nombre de cycles total à partir du nombre d'échantillons

        pour chaque cycle (1 à cyclesTotal inclus)
            on copie out dans in
            on met la largeur papillon comme étant cycle²
            pour chaque DTF à faire (de 0 à nSamples)
                Si le compteur de DTF est < la moitié de la largeur du papillon:
                    out[i] = DTF(in[i], in[i + cycle], compteurPap, largeurPap)
                Si le compteur de DTF est > que la moitié de la largeur du papillon:
                    out[i] = DTF(in[i - cycle], in[i], compteurPap, largeurPap)
                si compteurPap++ >= largeurPap
                    compteurPap = 0
    */
    unsigned int samples = data.size();
    vector<complex<double>> in(samples);
    for (int i = 0; i < data.size(); i++)
    {
        in[i] = data[mirror(i, 3)];
    }

    vector<complex<double>> out(samples);
    unsigned int maxExponent = 1;
    unsigned int butterflyOffset = 1;
    unsigned int nCycles = getNumberCycles(samples);
    for (int cycle = 1; cycle <= nCycles; cycle++)
    {
        if (cycle != 1)
            in = out;

        unsigned int butterflyWidth = pow(2.0, cycle);
        unsigned int butterflyCounter = 0;
        unsigned int exponent = 0;
        //cout << butterflyWidth << endl;
        for (int ftd = 0; ftd < samples; ftd++)
        {
            cout << "but " << butterflyCounter << endl;
            cout << "expo " << exponent << endl;
            //cout << (ftd > (butterflyWidth / 2)) << endl;
            if (butterflyCounter < (butterflyWidth / 2))
            {
                cout << "c1 " << ftd << " c2 " << (ftd + cycle) << endl;
                out[ftd] = iDFT(in[ftd], in[ftd + butterflyOffset], exponent, butterflyWidth, true);
                cout << "out" << out[ftd] << endl;
            }
            else if (butterflyCounter >= (butterflyWidth / 2))
            {
                cout << "c1 " << (ftd - cycle) << " c2 " << ftd << endl;
                out[ftd] = iDFT(in[ftd - butterflyOffset], in[ftd], exponent, butterflyWidth, false);
                cout << "out:" << out[ftd] << endl;
            }
            cout << endl;

            butterflyCounter++;
            if (butterflyCounter >= butterflyWidth)
            {
                butterflyCounter = 0;
            }

            exponent++;
            if (exponent >= maxExponent)
            {
                exponent = 0;
            }

        }
        cout << "--- CYCLE ---" << endl;

        maxExponent *= 2;

        butterflyOffset *= 2;
    }

    for (int i = 0; i < out.size(); i++)
    {
        out[i] = out[i] / (double)samples;
    }

    for (int i = 0; i < out.size(); i++)
    {
        cout << out[i] << endl;
    }
    cout << endl;
    for (int i = 0; i < in.size(); i++)
    {
        cout << in[i] << endl;
    }
    return out;
}

int main()
{
    vector<complex<double>> test1(8);
    test1[0] = complex<double>(1.0, 0.0);
    test1[1] = complex<double>(1.0, 0.0);
    test1[2] = complex<double>(1.0, 0.0);
    test1[3] = complex<double>(1.0, 0.0);
    test1[4] = complex<double>(1.0, 0.0);
    test1[5] = complex<double>(1.0, 0.0);
    test1[6] = complex<double>(1.0, 0.0);
    test1[7] = complex<double>(1.0, 0.0);

    vector<complex<double>> result1 = icooleyTurkey(test1);
    
}