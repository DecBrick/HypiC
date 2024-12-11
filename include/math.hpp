#pragma once
#include <vector>

namespace HypiC{
    double Forward_Difference(double f0, double f1, double f2, double x0, double x1, double x2);

    double Central_Difference(double f0, double f1, double f2, double x0, double x1, double x2);

    double Backward_Difference(double f0, double f1, double f2, double x0, double x1, double x2);

    //double Linear_Transition(double x, double cutoff, double L, double y1, double y2);

    std::vector<double> Thomas_Algorithm(std::vector<double> lower_diagonal, std::vector<double> diagonal, std::vector<double> upper_diagonal, std::vector<double> b);
}