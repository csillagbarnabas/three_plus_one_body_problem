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

template<typename T>
class orbit_vector{
	int N;
    std::vector<T> data;
	public:
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
    int Nsize()const{return N;}
    auto begin(){return data.begin();}
	auto cbegin() const{return data.cbegin();}
	auto end(){return data.end();}
	auto cend() const{return data.cend();}
};

template<typename T>
orbit_vector<T> operator+(orbit_vector<T> const& A, orbit_vector<T> const& B)
{
    orbit_vector<T> res{A.Nsize(),std::vector<T> (A.Nsize(), 0.0)};
    for(int i = 0; i<A.Nsize(); i++){res[i] = A[i]+B[i];}
	return res;
}

template<typename T>
orbit_vector<T> operator-(orbit_vector<T> const& A, orbit_vector<T> const& B)
{
    orbit_vector<T> res{A.Nsize(),std::vector<T> (A.Nsize(), 0.0)};
    for(int i = 0; i<A.Nsize(); i++){res[i] = A[i]-B[i];}
	return res;
}

template<typename T>
orbit_vector<T> operator*(orbit_vector<T> const& A,T const C)
{
    orbit_vector<T> res{A.Nsize(),std::vector<T> (A.Nsize(), 0.0)};
    for(int i = 0; i<A.Nsize(); i++){res[i] = A[i]*C;}
	return res;
}

template<typename T>
std::vector<T> acc(orbit_vector<T>& V){
    double x12 = V[1]-V[3];
    double y12 = V[2]-V[4];
    double rSquared_12 = x12*x12 + y12*y12;
    double rCubed_12 = rSquared_12 * std::sqrt(rSquared_12);
    double x13 = V[1]-V[5];
    double y13 = V[2]-V[6];
    double rSquared_13 = x13*x13 + y13*y13;
    double rCubed_13 = rSquared_13 * std::sqrt(rSquared_13);
    double x23 = V[3]-V[5];
    double y23 = V[4]-V[6];
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
    std::vector<T> acc{ax_1,ay_1,ax_2,ay_2,ax_3,ay_3};
    return acc;
}

template<typename T>
void accel(orbit_vector<T>& V,std::vector<T>& acce){
    double x12 = V[1]-V[3];
    double y12 = V[2]-V[4];
    double rSquared_12 = x12*x12 + y12*y12;
    double rCubed_12 = rSquared_12 * std::sqrt(rSquared_12);
    double x13 = V[1]-V[5];
    double y13 = V[2]-V[6];
    double rSquared_13 = x13*x13 + y13*y13;
    double rCubed_13 = rSquared_13 * std::sqrt(rSquared_13);
    double x23 = V[3]-V[5];
    double y23 = V[4]-V[6];
    double rSquared_23 = x23*x23 + y23*y23;
    double rCubed_23 = rSquared_23 * std::sqrt(rSquared_23);

    double ax12 = x12 / rCubed_12;
    double ax13 = x13 / rCubed_13;
    double ax23 = x23 / rCubed_23;

    double ay12 = y12 / rCubed_12;
    double ay13 = y13 / rCubed_13;
    double ay23 = y23 / rCubed_23;

    acce[0] = - G * (m2 * ax12 + m3 * ax13);
    acce[1] = - G * (m2 * ay12 + m3 * ay13);

    acce[2] = G * (m1 * ax12 - m3 * ax23);
    acce[3] = G * (m1 * ay12 - m3 * ay23);

    acce[4] = G * (m1 * ax13 + m2 * ax23);
    acce[5] = G * (m1 * ay13 + m2 * ay23);
}

template<typename T>
std::vector<T> acc4(orbit_vector<T>& V){
    double x14 = V[1]-V[13];
    double y14 = V[2]-V[14];
    double rSquared_14 = x14*x14 + y14*y14;
    double rCubed_14 = rSquared_14 * std::sqrt(rSquared_14);
    double x24 = V[3]-V[13];
    double y24 = V[4]-V[14];
    double rSquared_24 = x24*x24 + y24*y24;
    double rCubed_24 = rSquared_24 * std::sqrt(rSquared_24);
    double x34 = V[5]-V[13];
    double y34 = V[6]-V[14];
    double rSquared_34 = x34*x34 + y34*y34;
    double rCubed_34 = rSquared_34 * std::sqrt(rSquared_34);

    double ax14 = x14 / rCubed_14;
    double ax24 = x24 / rCubed_24;
    double ax34 = x34 / rCubed_34;

    double ay14 = y14 / rCubed_14;
    double ay24 = y24 / rCubed_24;
    double ay34 = y34 / rCubed_34;

    double ax4 = G * (m1 * ax14 + m2 * ax24 + m3 * ax34);
    double ay4 = G * (m1 * ay14 + m2 * ay24 + m3 * ay34);
    std::vector<T> acc4{ax4,ay4};
    return acc4;
}


template<typename T, typename F>
void velocity_verlet(orbit_vector<T>& A, F h,std::vector<T>& ac0,std::vector<T>& ac1){
    accel(A,ac0);
    std::vector<T> a41(2);
    if(A.Nsize() == 17){
        a41 = acc4(A);
        A[13] = A[13]+h*A[15]+h*h*a41[0]/2;
        A[14] = A[14]+h*A[16]+h*h*a41[1]/2;
    }
    for (int i = 0; i < 6; i++){A[i+1]  = A[i+1]+h*A[i+7]+h*h*ac0[i]/2;}
    accel(A,ac1);
    for (int i = 0; i < 6; i++){A[i+7] = A[i+7]+h*(ac0[i]+ac1[i])/2;}
    if(A.Nsize() == 17){
        std::vector<T> a42 = acc4(A);
        A[15] = A[15]+h*(a41[0] + a42[0])/2;
        A[16] = A[16]+h*(a41[1] + a42[1])/2;
    }
    A[0] = A[0] + h;
}

template<typename T>
void time_reversal(orbit_vector<T>& A){
    for (int i = 0; i < 6; i++){A[i+7] = - A[i+7];}
    A[15] = -A[15];
    A[16] = -A[16];
}

int main()
{
    auto t1 = std::chrono::high_resolution_clock::now();

    //my initial test parameters
    double v_add = 0.4;
    double x1 = 16, y1 = 5, x2 = -12, y2 = 0, x3 = 0, y3 = 20; // Astronomical unit - AU
    double v_x1=0.3+v_add,v_y1=1.02,v_x2=0.4+v_add,v_y2=1.1,v_x3=3.5+v_add,v_y3 = 0; // AU / Earth years
    double x4 = x2+1.5, y4 = y2, v_x4 = -0.55, v_y4 = -6;

    //Lagrange solution test
    //const double a = 20;
    //const double v0 = 0.5;
    //double x1=0,y1=a,x2=std::sqrt(3)*a/2,y2=-a/2,x3=-std::sqrt(3)*a/2,y3=-a/2; // Astronomical unit - AU
    //double v_x1=-v0,v_y1=0,v_x2=v0/2,v_y2=std::sqrt(3)*v0/2,v_x3=v0/2,v_y3=-std::sqrt(3)*v0/2; // AU / Earth years
    //double x4 = 0, y4 = y1 + 3, v_x4 = -6.5, v_y4 = 0;

    //orbit_vector<double> z{13,{0,x1,y1,x2,y2,x3,y3,v_x1,v_y1,v_x2,v_y2,v_x3,v_y3}};
    orbit_vector<double> z{17,{0,x1,y1,x2,y2,x3,y3,v_x1,v_y1,v_x2,v_y2,v_x3,v_y3,x4,y4,v_x4,v_y4}};
    double h = 1e-4;
    const std::string filename = "3p1body-VV2_myparams_tr_test.dat";
    //const std::string filename = "3body-RK4a_lagrange.dat";
    const int number_of_steps = 10000;
    std::vector<double> ac0(6);
    std::vector<double> ac1(6);
    std::cout << "3 body simulation\n";
    std::ofstream dataFile(filename);
    for (int step = 0; step < number_of_steps; step++){
        velocity_verlet(z,h,ac0,ac1);
        for (int j = 0; j < z.Nsize(); j++){
            dataFile << z[j] << '\t';}
        dataFile << '\n';}
    if(true){
        time_reversal(z);
        for (int step = 0; step < number_of_steps; step++){
            velocity_verlet(z,h,ac0,ac1);
            for (int j = 0; j < z.Nsize(); j++){
                dataFile << z[j] << '\t';}
            dataFile << '\n';}}
    dataFile.close();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Computation time:\n";
    std::cout << double((static_cast<std::chrono::duration<double, std::milli>>(t2-t1)).count()) <<" ms"<<std::endl;
}