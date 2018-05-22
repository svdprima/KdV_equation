#include <iostream>
#include <cstdlib>
#include "soliton.hpp"

int main (int argc, char** argv)
{
    if (argc == 1 || argc == 2)
    {
        std::cout << "Please, provide the number of points" << std::endl;
        exit (EXIT_FAILURE);
    }
    const double N = 2;
    const double T = 2;
    const unsigned int x_nu =  atoi(argv[1]);
    const unsigned int t_nu =  atoi(argv[2]);
    const double h = N / (double) (x_nu - 1); 
    const double tau = T / (double) (t_nu - 1);
    Korteweg_de_Vries S (x_nu, h, t_nu, tau, 0.022);
    S.Compute ();
    S.Print("output.txt");
    return 0;
}
