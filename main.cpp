#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>

const double G = 6.6743*1.98892*3.1536*3.1536/(1.496*1.496*1.496);
//gravitational constant in AU^3 / (Sun mass * Earth years^2)

const double m1 = 3, m2 = 2, m3 = 5; //in units of the Sun's mass

template<typename T, typename F>
std::vector<T> multiply(std::vector<T> & x, const F y){
    std::vector<T> res(x.size());
    for(int i = 0; i<x.size();i++){res[i] = x[i]*y;}
    return res;
}

template<typename T, typename F>
std::vector<T> add(std::vector<T> & x,std::vector<F> & y){
    std::vector<T> res(x.size());
    //for(int i = 0; i<x.size();i++){res[i] = x[i]+y[i];}
    std::transform(begin(x), end(x), begin(y),std::begin(res),[](T a,F b) { return a + b; });
    return res;
}

template<typename T>
class orbit_vector{
	int N;
    std::vector<T> data;
	public:
    //T *t, *x1, *y1, *v_x1, *v_y1, *x2, *y2, *v_x2, *v_y2, *x3, *y3, *v_x3, *v_y3;
    //t = &data[0], x1 = &data[1], y1 = &data[2], v_x1 = &data[3], v_y1 = &data[4];
    //x2 = &data[5], y2 = &data[6], v_x2 = &data[7], v_y2 = &data[8];
    //x3 = &data[9], y3 = &data[10], v_x3 = &data[11], v_y3 = &data[12];
    T&  operator()(int i)
    	{ return data[i]; }
    T const& operator()(int i) const
    	{ return data[i]; }
	T&  operator[](int i)
    { return data[i]; }
    T const& operator[](int i) const
    { return data[i]; }
    orbit_vector(): N(0), data(0) {};
	orbit_vector( orbit_vector const& ) = default;
	orbit_vector(orbit_vector&& m) : N{m.N}, data{std::move(m.data)} {m.N = 0;  };  
	orbit_vector<T>& operator=(orbit_vector const&) = default;
	orbit_vector<T>& operator=(orbit_vector && m){
		if(&data == &m.data){
			return *this;
		}
		else{
			N=m.N;
			data=std::move(m.data);
			m.N = 0; 
			return *this;
		}
	}
    orbit_vector( int n, std::vector<T> const& x ) : N(n), data(x){};
    int Nsize()const{
		return N;
		}

};

template<typename T>
orbit_vector<T> operator+(orbit_vector<T> const& A, orbit_vector<T> const& B)
{
    orbit_vector<T> res{A.Nsize(),std::vector<T> (A.Nsize(), 0.0)};
    for(int i = 0; i<A.Nsize(); i++){
        res[i] = A[i]+B[i];
    }
	return res;
}

template<typename T>
orbit_vector<T> operator*(orbit_vector<T> const& A,T const C)
{
    orbit_vector<T> res{A.Nsize(),std::vector<T> (A.Nsize(), 0.0)};
    for(int i = 0; i<A.Nsize(); i++){
        res[i] = A[i]*C;
    }
	return res;
}

template<typename T, typename F, typename RHS>
void RK4Step(std::vector<T>& z, F h, RHS derivs)
{
    std::vector<T> k1 = multiply(derivs(z),h);
    std::vector<T> k2 = multiply(derivs(add(z,multiply(k1,0.5))),h);
    std::vector<T> k3 = multiply(derivs(add(z,multiply(k2,0.5))),h);
    std::vector<T> k4 = multiply(derivs(add(z,k3)),h);
    z = add(z,multiply(add(add(k1, multiply(k2,2)), add(multiply(k3,2), k4)),1.0/6.0));
}

template<typename T, typename F, typename RHS>
void RK4Step_3body(orbit_vector<T>& z, F h, RHS derivs,
    orbit_vector<T>& k1,orbit_vector<T>& k2,orbit_vector<T>& k3,orbit_vector<T>& k4)
{
    k1 = derivs(z)*h;
    k2 = derivs((z+(k1*0.5)))*h;
    k3 = derivs((z+(k2*0.5)))*h;
    k4 = derivs((z+k3))*h;
    z = z+(k1 + (k2*2.0) + (k3*2.0)+ k4)*(1.0/6.0);
    //exit(-1);
}

template<typename T, typename F, typename L, typename RHS>
void adaptiveRK4Step(orbit_vector<T>& x, F tau, L accuracy, RHS derivs, orbit_vector<T>& x_half,
    orbit_vector<T>& x_full, orbit_vector<T>& Delta, orbit_vector<T>& scale,
    orbit_vector<T>& k1,orbit_vector<T>& k2,orbit_vector<T>& k3,orbit_vector<T>& k4)
{
    const double SAFETY = 0.9, PGROW = -0.2, PSHRINK = -0.25,
                 ERRCON = 1.89E-4, TINY = 1.0e-30;
    const int n = int(x.Nsize());
    for (int i = 0; i < n; i++)
        scale[i] = abs(x[i]) + abs(scale[i] * tau) + TINY;
    double err_max;
    while (true) {
        // take two half steps
        double tau_half = tau / 2;
        x_half = x;
        RK4Step_3body(x_half, tau_half, derivs,k1,k2,k3,k4);
        RK4Step_3body(x_half, tau_half, derivs,k1,k2,k3,k4);
        // take full step
        x_full = x;
        RK4Step_3body(x_full, tau, derivs,k1,k2,k3,k4);
        // estimate error
        for(int i = 0; i < n; i++){
            Delta[i] = x_half[i] - x_full[i];
        }
        err_max = 0;
        for (int i = 0; i < n; i++)
            err_max = std::max(err_max, std::abs(Delta[i]) / scale[i]);
        err_max /= accuracy;
        if (err_max <= 1.0)
            break;
        double tau_temp = SAFETY * tau * pow(err_max, PSHRINK);
        if (tau >= 0.0)
            tau = std::max(tau_temp, 0.1 * tau);
        else
            tau = std::min(tau_temp, 0.1 * tau);
        if (std::abs(tau) == 0.0) {
            std::cerr << "adaptiveRK4Step: step size underflow\naborting ..."
                 << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    tau *= (err_max > ERRCON ? SAFETY * pow(err_max, PGROW) : 5.0);
    for(int i=0; i < n; i++){
        x[i] = x_half[i] + Delta[i] / 15.0;
    }
}

//  Derivative vector for Newton's law of gravitation for 2 body:
std::vector<double> N2(const std::vector<double>& z) 
{
    double t = z[0], x1 = z[1], y1 = z[2], v_x1 = z[3], v_y1 = z[4];
    double x2 = z[5], y2 = z[6], v_x2 = z[7], v_y2 = z[8];
    double x12 = x1-x2;
    double y12 = y1-y2;
    double rSquared_ij = x12*x12 + y12*y12;
    double rCubed_ij = rSquared_ij * std::sqrt(rSquared_ij);
    double ax_ij = - G * m2 * x12 / rCubed_ij;
    double ay_ij = - G * m2 * y12 / rCubed_ij;
    std::vector<double> N2{1,v_x1,v_y1,ax_ij,ay_ij,v_x2,v_y2,- ax_ij*(m1/m2),- ay_ij*(m1/m2)};
    return N2;
}

//  Derivative vector for Newton's law of gravitation for 3 body:
std::vector<double> N3(const std::vector<double>& z) 
{
    double t = z[0], x1 = z[1], y1 = z[2], v_x1 = z[3], v_y1 = z[4];
    double x2 = z[5], y2 = z[6], v_x2 = z[7], v_y2 = z[8];
    double x3 = z[9], y3 = z[10], v_x3 = z[11], v_y3 = z[12];

    double x12 = x1-x2;
    double y12 = y1-y2;
    double rSquared_12 = x12*x12 + y12*y12;
    double rCubed_12 = rSquared_12 * std::sqrt(rSquared_12);
    double x13 = x1-x3;
    double y13 = y1-y3;
    double rSquared_13 = x13*x13 + y13*y13;
    double rCubed_13 = rSquared_13 * std::sqrt(rSquared_13);
    double x23 = x2-x3;
    double y23 = y2-y3;
    double rSquared_23 = x23*x23 + y23*y23;
    double rCubed_23 = rSquared_23 * std::sqrt(rSquared_23);

    double ax12 = x12 / rCubed_12;
    double ax13 = x13 / rCubed_13;
    double ax23 = x23 / rCubed_23;

    double ay12 = y12 / rCubed_12;
    double ay13 = y13 / rCubed_13;
    double ay23 = y23 / rCubed_23;

    double ax_1 = - G * (m2 * ax12 + m3 * ax13);
    double ay_1 = - G * (m2 * ay12 + m3 * ay13);

    double ax_2 = G * (m1 * ax12 - m3 * ax23);
    double ay_2 = G * (m1 * ay12 - m3 * ay23);

    double ax_3 = G * (m1 * ax13 + m2 * ax23);
    double ay_3 = G * (m1 * ay13 + m2 * ay23);
    
    std::vector<double> N3{1,v_x1,v_y1,ax_1,ay_1,v_x2,v_y2,ax_2,ay_2,v_x3,v_y3,ax_3,ay_3};
    return N3;
}

std::vector<double> N3_efficient(const std::vector<double>& z) 
{
    double x12 = z[1]-z[5];
    double y12 = z[2]-z[6];
    double rSquared_12 = x12*x12 + y12*y12;
    double rCubed_12 = rSquared_12 * std::sqrt(rSquared_12);
    double x13 = z[1]-z[9];
    double y13 = z[2]-z[10];
    double rSquared_13 = x13*x13 + y13*y13;
    double rCubed_13 = rSquared_13 * std::sqrt(rSquared_13);
    double x23 = z[5]-z[9];
    double y23 = z[6]-z[10];
    double rSquared_23 = x23*x23 + y23*y23;
    double rCubed_23 = rSquared_23 * std::sqrt(rSquared_23);

    double ax12 = x12 / rCubed_12;
    double ax13 = x13 / rCubed_13;
    double ax23 = x23 / rCubed_23;

    double ay12 = y12 / rCubed_12;
    double ay13 = y13 / rCubed_13;
    double ay23 = y23 / rCubed_23;

    double ax_1 = - G * (m2 * ax12 + m3 * ax13);
    double ay_1 = - G * (m2 * ay12 + m3 * ay13);

    double ax_2 = G * (m1 * ax12 - m3 * ax23);
    double ay_2 = G * (m1 * ay12 - m3 * ay23);

    double ax_3 = G * (m1 * ax13 + m2 * ax23);
    double ay_3 = G * (m1 * ay13 + m2 * ay23);
    
    std::vector<double> N3_efficient{1,z[3],z[4],ax_1,ay_1,z[7],z[8],ax_2,ay_2,z[11],z[12],ax_3,ay_3};
    return N3_efficient;
}

orbit_vector<double> N3o_efficient(const orbit_vector<double>& V) 
{
    double x12 = V[1]-V[5];
    double y12 = V[2]-V[6];
    double rSquared_12 = x12*x12 + y12*y12;
    double rCubed_12 = rSquared_12 * std::sqrt(rSquared_12);
    double x13 = V[1]-V[9];
    double y13 = V[2]-V[10];
    double rSquared_13 = x13*x13 + y13*y13;
    double rCubed_13 = rSquared_13 * std::sqrt(rSquared_13);
    double x23 = V[5]-V[9];
    double y23 = V[6]-V[10];
    double rSquared_23 = x23*x23 + y23*y23;
    double rCubed_23 = rSquared_23 * std::sqrt(rSquared_23);

    double ax12 = x12 / rCubed_12;
    double ax13 = x13 / rCubed_13;
    double ax23 = x23 / rCubed_23;

    double ay12 = y12 / rCubed_12;
    double ay13 = y13 / rCubed_13;
    double ay23 = y23 / rCubed_23;

    double ax_1 = - G * (m2 * ax12 + m3 * ax13);
    double ay_1 = - G * (m2 * ay12 + m3 * ay13);

    double ax_2 = G * (m1 * ax12 - m3 * ax23);
    double ay_2 = G * (m1 * ay12 - m3 * ay23);

    double ax_3 = G * (m1 * ax13 + m2 * ax23);
    double ay_3 = G * (m1 * ay13 + m2 * ay23);
    
    orbit_vector<double> N3o_efficient{V.Nsize(),{1,V[3],V[4],ax_1,ay_1,V[7],
    V[8],ax_2,ay_2,V[11],V[12],ax_3,ay_3}};
    return N3o_efficient;
}

int main()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    double x1 = 6, y1 = 0, x2 = -7, y2 = 0, x3 = 0, y3 = 20; // Astronomical unit - AU
    double v_x1=0.3,v_y1=1.02,v_x2=0.2,v_y2=-1.3,v_x3=4,v_y3 = 0; // AU / Earth years
    orbit_vector<double> z{13,{0,x1,y1,v_x1,v_y1,x2,y2,v_x2,v_y2,x3,y3,v_x3,v_y3}};
    const double h = 0.008;
    const double accuracy = 1e-7;
    const std::string filename = "3body-test.dat";
    const int number_of_steps = 10000;
    orbit_vector<double> k1{z.Nsize(),std::vector<double> (z.Nsize(), 0.0)};
    orbit_vector<double> k2{z.Nsize(),std::vector<double> (z.Nsize(), 0.0)};
    orbit_vector<double> k3{z.Nsize(),std::vector<double> (z.Nsize(), 0.0)};
    orbit_vector<double> k4{z.Nsize(),std::vector<double> (z.Nsize(), 0.0)};
    orbit_vector<double> x_half{z.Nsize(),std::vector<double> (z.Nsize(), 0.0)};
    orbit_vector<double> x_full{z.Nsize(),std::vector<double> (z.Nsize(), 0.0)};
    orbit_vector<double> Delta{z.Nsize(),std::vector<double> (z.Nsize(), 0.0)};
    orbit_vector<double> scale{z.Nsize(),std::vector<double> (z.Nsize(), 0.0)};
    std::cout << " 3+1 body simulation\n";
    std::ofstream dataFile(filename);
    for (int step = 0; step < number_of_steps; step++){
        //RK4Step_3body(z,h,N3o_efficient,k1,k2,k3,k4);
        adaptiveRK4Step(z, h, accuracy, N3o_efficient, x_half, x_full, Delta, scale, k1,k2,k3,k4);
        for (int j = 0; j < z.Nsize(); j++){
            dataFile << z[j] << '\t';}
        dataFile << '\n';}
    dataFile.close();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Computation time:\n";
    std::cout << double((static_cast<std::chrono::duration<double, std::milli>>(t2-t1)).count()) <<" ms"<<std::endl;
}