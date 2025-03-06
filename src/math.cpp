#include "math.hpp"
#include <iostream>
#include "HypiCpp.hpp"

namespace HypiC{

    double Backward_Difference(double f_l, double f_c, double x_l, double x_c){
        double h = x_c - x_l;
        return (f_c - f_l) / h;
        /*
        double h1 = x1 - x0;
        double h2 = x2 - x1;
        double c0 = -(2*h1+h2)/h1/(h1+h2);
        double c1 = (h1+h2)/(h1*h2);
        double c2 = -h1/h2/(h1+h2);
        return c0 * f0 + c1 * f1 + c2 * f2;*/
    }

    double Central_Difference(double f_l, double f_r, double x_l, double x_r){
        double h = x_r - x_l; //factor of 2 is included b/c we are using full cell centers 
        return (f_r - f_l) / h;
        /*
        double h1 = x1 - x0;
        double h2 = x2 - x1;
        double c0 = -h2/h1/(h1+h2);
        double c1 = -(h1+h2)/(h1*h2);
        double c2 = h1/h2/(h1+h2);
        return c0 * f0 + c1 * f1 + c2 * f2;*/
    }

    double Forward_Difference(double f_c, double f_r, double x_c, double x_r){
        double h = x_r - x_c;
        return (f_r - f_c) / h;
        /*        
        double h1 = x1 - x0;
        double h2 = x2 - x1;
        double c0 = h2/h1/(h1+h2);
        double c1 = -(h1+h2)/(h1*h2);
        double c2 = (h1+2*h2)/h2/(h1+h2);

        return c0 * f0 + c1 * f1 + c2 * f2;*/
    }


    /*
    double Linear_Transition(double x, double cutoff, double L, double y1, double y2){
        double x1 = cutoff - L/2;
        double x2 = cutoff + L/2;
        if (x < x1){
            return y1;
        } else if (x > x2){
            return y2;
        } else {
            double t = (x-x1)/(x2-x1);
            return t * (y2-y1) + y1;
        }
    }*/

    std::vector<double> Thomas_Algorithm(std::vector<double> lower_diagonal, std::vector<double> diagonal, std::vector<double> upper_diagonal, std::vector<double> b){
        //I pulled the description from Wiki, might need a better source
        size_t n;
        double w;
        //pull sizing information
        n = diagonal.size();
        
        //initialize result
        std::vector<double> x(n, 0.0);

        //#pragma omp parallel for private(w)
        //enter main loop 
        for (size_t i=1; i < n; i++){
            w = lower_diagonal[i-1] / diagonal[i-1];
            diagonal[i] -= w * upper_diagonal[i-1];
            b[i] -= w * b[i-1];
        }

        //back substitution 
        x[n-1] = b[n-1] / diagonal[n-1];
        //#pragma omp parallel for
        for (int i = n-2; i >= 0; --i){
            x[i] = (b[i] - upper_diagonal[i] * x[i+1]) / diagonal[i];
        }
        
        return x;
    }

}