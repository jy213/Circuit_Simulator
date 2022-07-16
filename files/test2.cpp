#include<iostream>
#include "part1.cpp"
#include<complex>
#include <string>
#include <eigen3/Eigen/Dense>

int main(){
    std::complex<double> t1 = (1,1);
    std::complex<double> t2 = (2,3);
    std::cout << t1 / t2;
}