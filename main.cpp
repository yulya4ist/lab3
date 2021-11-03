#include <iostream>
#include <fstream>
#include "cmath"

class table{
private:
    double* x;
    double* y;
    int n;
    double* a;
    double* b;
    double* c;
    double* d;
    bool use = false;
public:
    void refresh_spline(){
        if (use) {
            delete[] a;
            delete[] b;
            delete[] c;
            delete[] d;
        }
        use = true;
        a = new double[n];
        for (int i = 0; i < n; ++i) {
            a[i] = y[i];
        }
        b = new double[n-1];
        d = new double[n-1];
        auto* h = new double[n-1];
        for (int i = 0; i < n-1; ++i) {
            h[i] = x[i+1] - x[i];
        }
        auto* alpha = new double[n-1];
        for (int i = 1; i < n-1; ++i) {
            alpha[i] = (3/h[i])*(a[i+1]-a[i])-(3/h[i-1])*(a[i]-a[i-1]);
        }
        c = new double[n];
        auto* l = new double[n];
        auto* nu = new double[n];
        auto* z = new double[n];
        l[0] = 1; nu[0] = 0; z[0] = 0;
        for (int i = 0; i < n-1; ++i) {
            l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*nu[i-1];
            nu[i] = h[i]/l[i];
            z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
        }
        l[n-1] = 1; z[n-1] = 0; c[n-1] = 0;
        for (int i = n-2; i >= 0; --i) {
            c[i] = z[i]-nu[i]*c[i+1];
            b[i] = (a[i+1]-a[i])/h[i] - h[i]*(c[i+1]+2*c[i])/3;
            d[i] = (c[i+1]-c[i])/(3*h[i]);
        }
        delete[] alpha;
        delete[] l;
        delete[] nu;
        delete[] z;
    }
    table() {
        x = new double[1];
        y = new double[1];
        x[0] = 0;
        y[0] = 0;
        n = 1;
        this->refresh_spline();
    }
    table(double (*fun1)(double), double a, double b, int n_add) {
        double step = (b-a)/n_add;
        double x_curr = a;
        x = new double[n_add+1];
        y = new double[n_add+1];
        for (int i = 0; i <= n_add; ++i) {
            x[i] = x_curr;
            y[i] = fun1(x_curr);
            x_curr += step;
        }
        n = n_add+1;
        this->refresh_spline();
    }
    table(double (*fun1)(double), double a, double b, int n_add, bool Cheb) {
        double x_curr;
        x = new double[n_add+1];
        y = new double[n_add+1];
        for (int i = 0; i <= n_add; ++i) {
            x_curr = ((a+b)/2)+((b-a)/2)*std::cos((2*i+1)*M_PI/(2*(n_add+1)));
            x[i] = x_curr;
            y[i] = fun1(x_curr);
        }
        n = n_add+1;
        this->refresh_spline();
    }
    table(table const &A, double a, double b, int n_add, bool lagrange) {
        double step = (b-a)/n_add;
        double x_curr = a;
        x = new double[n_add+1];
        y = new double[n_add+1];
        for (int i = 0; i <= n_add; ++i) {
            x[i] = x_curr;
            if (lagrange) {
                y[i] = A.LagrangeInterp(x_curr);
            } else
            {
                y[i] = A.SplainInterp(x_curr);
            }
            x_curr += step;
        }
        n = n_add+1;
        this->refresh_spline();
    }
    void push_back(double x_add, double y_add) {
        auto* x_new = new double[n+1];
        auto* y_new = new double[n+1];
        for (int i = 0; i < n; ++i) {
            x_new[i] = x[i];
            y_new[i] = y[i];
        }
        x_new[n] = x_add;
        y_new[n] = y_add;
        double* c;
        c = x_new;
        x_new = x;
        x = c;
        delete[] x_new;
        c = y_new;
        y_new = y;
        y = c;
        delete[] y_new;
        n++;
        this->refresh_spline();
    }
    void output_table(){
        std::cout << std::endl;
        std::cout << "N  |";
        for (int i = 0; i < n; ++i) {
            std::cout << "| ";
            std::cout.width(13);
            std::cout<< i;
            std::cout << " |";
        }
        std::cout << "|" << std::endl;
        std::cout << "x  |";
        for (int i = 0; i < n; ++i) {
            std::cout << "| ";
            std::cout.width(13);
            std::cout<< x[i];
            std::cout << " |";
        }
        std::cout << "|" << std::endl;
        std::cout << "y  |";
        for (int i = 0; i < n; ++i) {
            std::cout << "| ";
            std::cout.width(13);
            std::cout<< y[i];
            std::cout << " |";
        }
        std::cout << "|" << std::endl;
    }
    [[nodiscard]] double LagrangeInterp(double x_in) const {
        double result = 0;
        for (int i = 0; i < n; ++i) {
            double ck = 1;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    ck = ck*(x_in-x[j])/(x[i]-x[j]);
                }
            }
            result = result + ck*y[i];
        }
        return result;
    }
    [[nodiscard]] double SplainInterp(double x_in) const {
        int i = 1;
        while ((x_in > x[i])&&(i<n-1)) {
            i++;
        }
        return a[i-1]+b[i-1]*(x_in-x[i-1])+c[i-1]*std::pow(x_in-x[i-1],2)+d[i-1]*std::pow(x_in-x[i-1],3);
    }
    ~table() {
        delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
        delete[] x;
        delete[] y;
    }
};


double var10(double x) {
    double r1 = std::sin((std::pow(x,4)+std::pow(x,3)-3*x+3-std::pow(30,1/3))/2);
    double r2 = std::tanh((4*std::sqrt(3)*std::pow(x,3)-2*x-6*std::sqrt(2)+1)/
            (-2*std::sqrt(3)*std::pow(x,3)+x+3*std::sqrt(2)));
    return r1+r2+1.2;
}

double var21(double x) {
    double r1 = std::pow((std::sqrt(3)*std::pow(x,3)-2*x+5)/(7+std::sqrt(7)),std::log(3)/std::log(10));
    double r2 = std::asin((std::pow(x,2)+x+std::sqrt(3))/(2*x-2));
    return r1+r2;
}

double var13(double x) {
    double r1 = std::sin((2*std::pow(x,2)-x+2*std::pow(7,1/3)-5)/2);
    double r2 = std::exp((std::pow(x,2)+2*x+1)/(7*x+1));
    return r1+r2-1.5;
}

double var0(double x) {
    return 1;
}

void lab3(double (*fun1)(double), int a, int b, std::string name){
    std::ofstream out1;
    out1.open("/home/qw/Рабочий стол/laba3_2/"+name+".txt");
    std::ofstream out_x1;
    out_x1.open("/home/qw/Рабочий стол/laba3_2/x_1"+name+".txt");
    std::ofstream out_y1;
    out_y1.open("/home/qw/Рабочий стол/laba3_2/y_1"+name+".txt");
    std::ofstream out_y2;
    out_y2.open("/home/qw/Рабочий стол/laba3_2/y_2"+name+".txt");
    table T(fun1,a,b,4);
    double err = 0;
    for (double i = a; i < b; i=i+0.001) {
        double f2 = T.LagrangeInterp(i);
        double f1 = fun1(i);
        if (std::abs(f2-f1)>err) err = std::abs(f2-f1);
        out_x1 << i << std::endl;
        out_y1 << f1 << std::endl;
        out_y2 << f2 << std::endl;
    }
    out_x1 << b << std::endl;
    out_y1 << T.LagrangeInterp(b);
    out_y2 << fun1(b);
    out1<< 4 << "   " << err;
    table T12(fun1,a,b,4, true);
    err = 0;
    for (double i = a; i < b; i=i+0.001) {
        double f1 = T.LagrangeInterp(i);
        double f2 = fun1(i);
        if (std::abs(f2-f1)>err) err = std::abs(f2-f1);
    }
    out1 << "   " << err << std::endl;
    for (int i = 8; i < 129; i=i*2) {
        table T1(fun1,a,b,i);
        table T2(fun1,a,b,i,true);
        double err1 = 0;
        double err2 = 0;
        for (double j = a; j < b; j=j+0.001) {
            double f1 = T1.LagrangeInterp(j);
            double f2 = fun1(j);
            double f3 = T2.LagrangeInterp(j);
            if (std::abs(f2-f1)>err1) err1 = std::abs(f2-f1);
            if (std::abs(f2-f3)>err2) err2 = std::abs(f2-f3);
        }
        out1 << i << "   " << err1 << "   " << err2 << std::endl;
    }
    table T_many(fun1,a,b,128);
    std::ofstream out_x1_128;
    out_x1_128.open("/home/qw/Рабочий стол/laba3_2/x_1_128"+name+".txt");
    std::ofstream out_y1_128;
    out_y1_128.open("/home/qw/Рабочий стол/laba3_2/y_1_128"+name+".txt");
    std::ofstream out_y2_128;
    out_y2_128.open("/home/qw/Рабочий стол/laba3_2/y_2_128"+name+".txt");
    for (double i = a; i < b; i=i+0.001) {
        out_x1_128 << i << std::endl;
        out_y1_128 << fun1(i) << std::endl;
        out_y2_128 << T_many.LagrangeInterp(i) << std::endl;
    }
    out_x1_128 << b;
    out_y1_128 << fun1(b);
    out_y2_128 << T_many.LagrangeInterp(b);
    std::ofstream out_x2;
    out_x2.open("/home/qw/Рабочий стол/laba3_2/x_2"+name+".txt");
    std::ofstream out_y2_1;
    out_y2_1.open("/home/qw/Рабочий стол/laba3_2/y_2_1"+name+".txt");
    std::ofstream out_y2_2;
    out_y2_2.open("/home/qw/Рабочий стол/laba3_2/y_2_2"+name+".txt");
    table TCheb1(fun1,a,b,4,true);
    for (double i = a; i < b; i=i+0.001) {
        out_x2 << i << std::endl;
        out_y2_1 << fun1(i) << std::endl;
        out_y2_2 << TCheb1.LagrangeInterp(i) << std::endl;
    }
    out_x2 << b;
    out_y2_1 << fun1(b);
    out_y2_2 << TCheb1.LagrangeInterp(b);
    std::ofstream out_x2_128;
    out_x2_128.open("/home/qw/Рабочий стол/laba3_2/x_2_128"+name+".txt");
    std::ofstream out_y2_1_128;
    out_y2_1_128.open("/home/qw/Рабочий стол/laba3_2/y_2_1_128"+name+".txt");
    std::ofstream out_y2_2_128;
    out_y2_2_128.open("/home/qw/Рабочий стол/laba3_2/y_2_2_128"+name+".txt");
    table TCheb1_128(fun1,a,b,128,true);
    for (double i = a; i < b; i=i+0.001) {
        out_x2_128 << i << std::endl;
        out_y2_1_128 << fun1(i) << std::endl;
        out_y2_2_128 << TCheb1_128.LagrangeInterp(i) << std::endl;
    }
    out_x2_128 << b;
    out_y2_1_128 << fun1(b);
    out_y2_2_128 << TCheb1_128.LagrangeInterp(b);
    std::ofstream out_x1_splain;
    out_x1_splain.open("/home/qw/Рабочий стол/laba3_2/x_1_splain"+name+".txt");
    std::ofstream out_y1_splain;
    out_y1_splain.open("/home/qw/Рабочий стол/laba3_2/y_1_splain"+name+".txt");
    std::ofstream out_y2_splain;
    out_y2_splain.open("/home/qw/Рабочий стол/laba3_2/y_2_splain"+name+".txt");
    for (double i = a; i < b; i=i+0.001) {
        out_x1_splain << i << std::endl;
        out_y1_splain << fun1(i) << std::endl;
        out_y2_splain << T.SplainInterp(i) << std::endl;
    }
    out_x1_splain << b;
    out_y1_splain << fun1(b);
    out_y2_splain << T.SplainInterp(b);
    std::ofstream out_x1_splain_64;
    out_x1_splain_64.open("/home/qw/Рабочий стол/laba3_2/x_1_splain_64"+name+".txt");
    std::ofstream out_y1_splain_64;
    out_y1_splain_64.open("/home/qw/Рабочий стол/laba3_2/y_1_splain_64"+name+".txt");
    std::ofstream out_y2_splain_64;
    out_y2_splain_64.open("/home/qw/Рабочий стол/laba3_2/y_2_splain_64"+name+".txt");
    for (double i = a; i < b; i=i+0.001) {
        out_x1_splain_64 << i << std::endl;
        out_y1_splain_64 << fun1(i) << std::endl;
        out_y2_splain_64 << T_many.SplainInterp(i) << std::endl;
    }
    out_x1_splain_64 << b;
    out_y1_splain_64 << fun1(b);
    out_y2_splain_64 << T_many.SplainInterp(b);
};

int main() {
    lab3(var13,0,1,"Kiseleva");
    //lab3(var0,-1,1,"Test");
    return 0;
}
