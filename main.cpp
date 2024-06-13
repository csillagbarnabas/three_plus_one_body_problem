#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>

const double G = 6.6743*1.98892*3.1536*3.1536/(1.496*1.496*1.496);
//gravitational constant in AU^3 / (Sun mass * Earth years^2)

template<typename T, typename F>
std::vector<T> multiply(std::vector<T> & x, const F y){
    std::vector<T> res(x.size());
    for(int i = 0; i<x.size();i++){
        res[i] = x[i]*y;
    }
    return res;
}

template<typename T, typename F>
std::vector<T> add(std::vector<T> & x,std::vector<F> & y){
    std::vector<T> res(x.size());
    for(int i = 0; i<x.size();i++){
        res[i] = x[i]+y[i];
    }
    return res;
}

template<typename T, typename F, typename RHS>
void RK4Step(std::vector<T>& z, F h, RHS derivs,const double m1,const double m2)
{
    std::vector<T> k1 = multiply(derivs(z,m1,m2),h);
    std::vector<T> k2 = multiply(derivs(add(z,multiply(k1,0.5)),m1,m2),h);
    std::vector<T> k3 = multiply(derivs(add(z,multiply(k2,0.5)),m1,m2),h);
    std::vector<T> k4 = multiply(derivs(add(z,k3),m1,m2),h);
    z = add(z,multiply(add(add(k1, multiply(k2,2)), add(multiply(k3,2), k4)),1.0/6.0));
}

//  Derivative vector for Newton's law of gravitation for 2 body:
std::vector<double> N2(const std::vector<double>& z,const double m1,const double m2) 
{
    double t = z[0], x1 = z[1], y1 = z[2], v_x1 = z[3], v_y1 = z[4];
    double x2 = z[5], y2 = z[6], v_x2 = z[7], v_y2 = z[8];
    double x12 = x1-x2;
    double y12 = y1-y2;
    double rSquared_ij = x12*x12 + y12*y12;
    double rCubed_ij = rSquared_ij * std::sqrt(rSquared_ij);
    double Fx_ij = - G * m2 * x12 / rCubed_ij;
    double Fy_ij = - G * m2 * y12 / rCubed_ij;
    std::vector<double> N2{1,v_x1,v_y1,Fx_ij,Fy_ij,v_x2,v_y2,- Fx_ij*(m1/m2),- Fy_ij*(m1/m2)};
    return N2;
}

int main()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    double x1 = 6, y1 = 0, x2 = -7, y2 = 0; // Astronomical unit - AU
    double v_x1 = 0.3, v_y1 = 1.02, v_x2 = 0.2, v_y2 = -1.3; // AU / Earth years
	std::vector<double> z{0,x1,y1,v_x1,v_y1,x2,y2,v_x2,v_y2};
    const double h = 0.005;
    const double m1 = 3;
    const double m2 = 2;
    const std::string filename = "2body-test.dat";
    const int number_of_steps = 10000;

    std::cout << " 3+1 body simulations\n";
    std::ofstream dataFile(filename);
    for (int step = 0; step < number_of_steps; step++){
        RK4Step(z,h,N2,m1,m2);
        for (int i = 0; i < z.size(); i++){
            dataFile << z[i] << '\t';}
        dataFile << '\n';}
    dataFile.close();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Computation time:\n";
    std::cout << double((static_cast<std::chrono::duration<double, std::milli>>(t2-t1)).count()) <<" ms"<<std::endl;
}