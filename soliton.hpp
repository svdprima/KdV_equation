#ifndef SOLITON_HPP
#define SOLITON_HPP

#include <vector>
#include <cmath>
#include <fstream>
#include <time.h>
#include <cstdlib>

class Korteweg_de_Vries
{
private:
    double x_nu;
    double t_nu;
    double delta;
    double h;
    double tau;
    std::vector <double> u;
public:
    Korteweg_de_Vries (double x_pts, double x_step, double t_pts, double t_step, double d)
    {
        delta = d;
        x_nu = x_pts;
        if (x_nu <= 5)
        {
            std::cout << "Not enough points in spacial domain" << std::endl;
            exit (-1);
        }
        t_nu = t_pts;
        if (t_nu <= 3)
        {
            std::cout << "Not enough points in time domain" << std::endl;
        }
        h = x_step;
        tau = t_step;
        if (tau / pow(h, 3) > 1 / (4 + 2 * pow(u0, 2) * pow(h, 2)))
        {
            std::cout << "The method is unstable, try other parameters" << std::endl;
            exit (EXIT_SUCCESS);
        }
        if (tau >= 1 / M_PI)
        {
            std::cout << "Try more points on time domain" << std::endl;
            exit (EXIT_SUCCESS);
        }
        u = std::vector <double> (x_nu * t_nu, 0);
    }
    void Compute ()
    {
        for (unsigned int i = 0; i < x_nu; i++)
        {
            u[i]        = cos(M_PI * i * h);
            u[i + x_nu] = cos(M_PI * i * h) / (1 - sin(M_PI * i * h) * M_PI * tau); 
        }
        for (unsigned int j = 1; j < t_nu - 1; j++)
            for (unsigned int i = 0; i < x_nu; i++)
            {
                if (i == 0)
                {
                    u [i + x_nu * (j + 1)] = u [i + x_nu * (j - 1)] - tau / 3 / h * (u [i + 1 + x_nu * j] + u [i + x_nu * j] + u[x_nu - 1 + x_nu * j]) *
                                             (u [i + 1 + x_nu * j] - u [x_nu - 1 + x_nu * j]) - pow(delta, 2) * tau / pow(h, 3) * 
                                             (u [i + 2 + x_nu * j] - 2 * u [i + 1 + x_nu * j] + 2 * u [x_nu - 1 + x_nu * j] - u [x_nu - 2 + x_nu * j]);
                }
                if (i == 1)
                {
                    u [i + x_nu * (j + 1)] = u [i + x_nu * (j - 1)] - tau / 3 / h * (u [i + 1 + x_nu * j] + u [i + x_nu * j] + u[i - 1 + x_nu * j]) *
                                             (u [i + 1 + x_nu * j] - u [i - 1 + x_nu * j]) - pow(delta, 2) * tau / pow(h, 3) * 
                                             (u [i + 2 + x_nu * j] - 2 * u [i + 1 + x_nu * j] + 2 * u [i - 1 + x_nu * j] - u [x_nu - 1 + x_nu * j]);
                }
                if (i == x_nu - 1)
                {
                    u [i + x_nu * (j + 1)] = u [i + x_nu * (j - 1)] - tau / 3 / h * (u [0 + x_nu * j] + u [i + x_nu * j] + u[i - 1 + x_nu * j]) *
                                             (u [0 + x_nu * j] - u [i - 1 + x_nu * j]) - pow(delta, 2) * tau / pow(h, 3) * 
                                             (u [1 + x_nu * j] - 2 * u [0 + x_nu * j] + 2 * u [i - 1 + x_nu * j] - u [i - 2 + x_nu * j]);
                }
                if (i == x_nu - 2)
                {
                    u [i + x_nu * (j + 1)] = u [i + x_nu * (j - 1)] - tau / 3 / h * (u [i + 1 + x_nu * j] + u [i + x_nu * j] + u[i - 1 + x_nu * j]) *
                                             (u [i + 1 + x_nu * j] - u [i - 1 + x_nu * j]) - pow(delta, 2) * tau / pow(h, 3) * 
                                             (u [0 + x_nu * j] - 2 * u [i + 1 + x_nu * j] + 2 * u [i - 1 + x_nu * j] - u [i - 2 + x_nu * j]);
                }
                if (i != 0 && i != 1 && i != x_nu - 1 && i != x_nu - 2)
                {
                    u [i + x_nu * (j + 1)] = u [i + x_nu * (j - 1)] - tau / 3 / h * (u [i + 1 + x_nu * j] + u [i + x_nu * j] + u[i - 1 + x_nu * j]) *
                                             (u [i + 1 + x_nu * j] - u [i - 1 + x_nu * j]) - pow(delta, 2) * tau / pow(h, 3) * 
                                             (u [i + 2 + x_nu * j] - 2 * u [i + 1 + x_nu * j] + 2 * u [i - 1 + x_nu * j] - u [i - 2 + x_nu * j]);
                }
                if (u[i + x_nu * (j + 1)] == INFINITY || u[i + x_nu * (j + 1)] == -INFINITY)
                {
                    std::cout << "The method started to diverge" << std::endl;
                    exit (-1);
                }
            }
    }
    void Print (const char* filename)
    {
        std::ofstream output;
        output.precision (20);
        output.open (filename);
        for (unsigned int i = 0; i < t_nu; i+= 1000) //since the number of points on t-axis is extremly high, 
            for (unsigned int j = 0; j < x_nu; j++)  //the function will create extremly huge files with usual i++
                output << std::fixed << j * h << " " << u[j + x_nu * i] << " " << i * tau << std::endl;
        output.close ();
    }
};

#endif
