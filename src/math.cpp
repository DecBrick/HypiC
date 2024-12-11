#include "math.hpp"

namespace HypiC{

    double Forward_Difference(double f0, double f1, double f2, double x0, double x1, double x2){
        double h1 = x1 - x0;
        double h2 = x2 - x1;
        double c0 = -(2*h1+h2)/h1/(h1+h2);
        double c1 = (h1+h2)/(h1*h2);
        double c2 = -h1/h2/(h1+h2);
        return c0 * f0 + c1 * f1 + c2 * f2;
    }

    double Central_Difference(double f0, double f1, double f2, double x0, double x1, double x2){
        double h1 = x1 - x0;
        double h2 = x2 - x1;
        double c0 = -h2/h1/(h1+h2);
        double c1 = -(h1+h2)/(h1*h2);
        double c2 = h1/h2/(h1+h2);
        return c0 * f0 + c1 * f1 + c2 * f2;
    }

    double Backward_Difference(double f0, double f1, double f2, double x0, double x1, double x2){
        double h1 = x1 - x0;
        double h2 = x2 - x1;
        double c0 = h2/h1/(h1+h2);
        double c1 = -(h1+h2)/(h1*h2);
        double c2 = (h1+2*h2)/h2/(h1+h2);

        return c0 * f0 + c1 * f1 + c2 * f2;
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

        //enter main loop 
        for (size_t i=1; i < n; i++){
            w = lower_diagonal[i-1] / diagonal[i-1];
            diagonal[i] -= w * upper_diagonal[i-1];
            b[i] -= w * b[i-1];
        }

        //back substitution 
        x[n-1] = b[n-1] / diagonal[n-1];
        int jj = n-2;
        for (size_t i = n-2; i-- > 0; ){
            x[i] = (b[i] - upper_diagonal[i] * x[i+1]) / diagonal[i];
        }
        
        return x;
    }

}