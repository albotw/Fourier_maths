#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <complex>

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
unsigned long mirror(unsigned long data, unsigned char nBits)
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

    return i+1;
}
/*Mini DFT entre deux valeurs*/
/*Tableau de 2 complexes, échantillon actuel, nombre d'échantillons*/
complex<double> DFT(complex<double> c1, complex<double> c2, int k, int n)
{
    complex<double> omega = exp(-(2 * M_PI * I * (double)k) / (double)n);
    return c1 + omega * c2;
}

void cooleyTurkey()
{
    /*
        On copie image dans out avec les inversions de bit.
        On calcule le nombre de cycles total à partir du nombre d'échantillons
        pour chaque cycle (1 à cyclesTotal inclus)
            on met la largeur papillon comme étant cycle²
            pour chaque DTF à faire (de 0 à nSamples)
                Si le compteur de DTF est < la moitié de la largeur du papillon:
                    out[i] = DTF(out[i], out[i + cycle], 0, largeurPap)
                Si le compteur de DTF est > que la moitié de la largeur du papillon:
                    out[i] = DTF(out[i - cycle], out[i], 1, largeurPap) 
    */
}

int main()
{
    printf("Hello, World !\n");
    unsigned char nombre = 0b101;
    affBin(nombre, 3, '0', '1');
    affBin(mirror(nombre, 3), 3, '0', '1');
}