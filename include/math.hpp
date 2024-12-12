#pragma once
#include <vector>

namespace HypiC{
    double Forward_Difference(double f_c, double f_r, double x_c, double x_r);

    double Central_Difference(double f_l, double f_r, double x_l, double x_r);

    double Backward_Difference(double f_l, double f_c, double x_l, double x_c);

    //double Linear_Transition(double x, double cutoff, double L, double y1, double y2);

    std::vector<double> Thomas_Algorithm(std::vector<double> lower_diagonal, std::vector<double> diagonal, std::vector<double> upper_diagonal, std::vector<double> b);
}