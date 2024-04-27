/*#include <iostream>
#include <cmath>*/


//1
/*int main()
{
	using namespace std;
	
	double x, y, xn, xk, dx, i;
	
	cout << "Vvedite xn \n";
	cin >> xn;
	cout << "Vvedite xk \n";
	cin >> xk;
	cout << "Vvedite dx \n";
	cin >> dx;
	double n = (xk - xn) / dx + 1;
	double s = 0;
	double p = 1;
	for (i = 0; i < n; ++i)
	{
		x = xn + i * dx;
		y = x * x * exp(sin(x) + cos(x));
		s += y;
		p *= y;
		cout <<"x = "<<x<<" y = "<<y<<"\n";
	}
	cout << "sum: " <<s<< " geometr: "<< p<<'\n';
	s /= n;
	p = pow(p,double(1/n));
	cout<<"Avarage sum:"<<s<< endl;
	cout<<"Avarage:"<<p<< endl;

}*/


//2
/*using namespace std;

int main() {
	double xn, xk, dx, x, y, ymax, ymin, xmax, xmin;
	cout << "Vvedite xn: ";
	cin >> xn;
	cout << "Vvedite xk: ";
	cin >> xk;
	cout << "Vvedite dx: ";
	cin >> dx;
	ymax = ymin = xn * sin(xn) - xn / 2 * cos(xn / 2);
	xmax = xmin = xn;
	for (x = xn; x <= xk; x += dx) {
		y = x * sin(x) - x / 2 * cos(x / 2);
		if (y > ymax) {
			ymax = y;
			xmax = x;
		}
		if (y < ymin) {
			ymin = y;
			xmin = x;

			cout << "Y: " << y << ", X = " << x << endl;
		}
	}
	cout << "Ymax: " << ymax << ", Xmax = " << xmax << endl;
	cout << "Ymin: " << ymin << ", Xmin = " << xmin << endl;
	return 0;
}*/


//3

/*using namespace std;

int main()
{
    double x, xn, xk, dx, f1, f2, f3, min, x1, x2, xe1, xe2, dxe;
    int i, n, k;

    cout << "Enter xn, xk, dx: ";
    cin >> xn >> xk >> dx;

    n = (xk - xn) / dx + 1;
    x = xn;
    f2 = x * sin(3 * x) * cos(x);
    cout << "x=" << x << " f=" << f2 << endl;

    x += dx;
    f3 = x * sin(3 * x) * cos(x);
    cout << "x=" << x << " f=" << f3 << endl;

    k = 0;
    min = xk - xn;

    for (i = 2; i < n; i++)
    {
        f1 = f2;
        f2 = f3;
        x = xn + i * dx;
        f3 = x * sin(3 * x) * cos(x);
        cout << "x=" << x << " f=" << f3 << endl;

        if (f2 > f1 && f2 > f3)
        {
            k++;
            cout << "k=" << k << " x=" << x - dx << " f=" << f2 << " max" << endl;

            if (k == 1)
            {
                x1 = x - dx;
            }

            if (k > 1)
            {
                x2 = x - dx;
                dxe = x2 - x1;

                if (dxe < min)
                {
                    min = dxe;
                    xe1 = x1;
                    xe2 = x2;
                }

                x1 = x2;
            }
        }

        if (f2 < f1 && f2 < f3)
        {
            k++;
            cout << "k=" << k << " x=" << x - dx << " f=" << f2 << " min" << endl;

            if (k == 1)
            {
                x1 = x - dx;
            }

            if (k > 1)
            {
                x2 = x - dx;
                dxe = x2 - x1;

                if (dxe < min)
                {
                    min = dxe;
                    xe1 = x1;
                    xe2 = x2;
                }

                x1 = x2;
            }
        }
    }

    cout << "min=" << min << " xe1=" << xe1 << " xe2=" << xe2 << endl;
    return 0;
}
*/


//dz1
/*
using namespace std;

int main() {
    const int n = 25;
    int X[n];

    
    for (int i = 0; i < n; i++) {
        X[i] = rand() % 101 - 50;
    }

    
    cout << "Dannyi massiv: ";
    for (int i = 0; i < n; i++) {
        cout << X[i] << " ";
    }
    cout << endl;

    int first_positive = -1; 
    int last_positive = -1; 
    for (int i = 0; i < n; i++) {
        if (X[i] > 0) {
            if (first_positive == -1) {
                first_positive = i;
            }
            last_positive = i;
        }
    }

    
    if (first_positive == -1) {
        cout << "V massive net poloshitelnogo elementa" << endl;
        return 0;
    }

    
    cout << "Pervyi poloshitelniu element: " << X[first_positive] << endl;
    cout << "Posledniu poloshitelniu element: " << X[last_positive] << endl;

    int sum = 0; 
    for (int i = first_positive + 1; i < last_positive; i++) {
        sum += X[i];
    }
    cout << "Summa chisel  meshdy pervyv i poslednim poloshitelnym elementom : " << sum << endl;

    int product = 1; 
    for (int i = 0; i < first_positive; i++) {
        product *= X[i];
    }
    for (int i = last_positive + 1; i < n; i++) {
        product *= X[i];
    }
    cout << "Proizvedenie : " << product << endl;

    return 0;
}*/


//dz2
/*
using namespace std;

double func(double x) {
    return x * (sin(3 * x) + cos(2 * x));
}
int main() {
    double x = -3.5;
    double step = 0.001;
    double max_dist = 0.0;
    double prev_extrimum =func(x);
    double curr_extrimum, curr_dist;
    bool increasing = false;


    while (x <= 4.5) {
        double y = func(x);
        
        if (y > prev_extrimum && increasing) {
            curr_extrimum = prev_extrimum;
            curr_dist = y - prev_extrimum;
            if (curr_dist > max_dist) {
                max_dist = curr_dist;
            }
            increasing = false;
        }
        else if (y < prev_extrimum && !increasing) {
            curr_extrimum = prev_extrimum;
            curr_dist = prev_extrimum - y;
            if (curr_dist > max_dist) {
                max_dist = curr_dist;
            }
            increasing = true;
        }
        prev_extrimum = y;
        x += step;
    }
    cout << "Max distance " << max_dist << endl;
}*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
#include <iostream>
#include <vector>
#include <array>
#include <cassert>
#include <cmath>

double f(double x)
{
    return x + std::log(x);
}

void dichotomy()
{
    double a = 0.1;
    double b = 1.0;

    std::vector<double> epsRange = { 0.000001, 0.0000001, 0.00000001 };

    std::cout << "Going dichotomy:" << std::endl;

    for (auto eps : epsRange)
    {
        int m = 0;
        double L;

        do
        {
            double xm = (a + b) / 2.0;
            L = (b - a);

            double fxm = f(xm);

            double x1 = a + L / 4.0;
            double x2 = b - L / 4.0;

            double fx1 = f(x1);
            double fx2 = f(x2);

            if (fx1 < fxm)
            {
                b = xm;
                xm = x1;
            }
            else
            {
                if (fx2 < fxm)
                {
                    a = xm;
                    xm = x2;
                }
                else
                {
                    a = x1;
                    b = x2;
                }
            }

            m++;
        } while (L > eps);

        std::cout << "For eps = " << eps << " m = " << m << std::endl;
    }
}

int main()
{
    dichotomy();
    return 0;
}*/


/*
#include <iostream>
#include <vector>
#include <array>
#include <cassert>
#include <cmath>

double f(double x)
{
    return x + std::log(x);
}

void goldenRatio()
{
    double a = 0.1;
    double b = 1.0;

    std::vector<double> epsRange = { 0.000001, 0.0000001, 0.00000001 };

    std::cout << "Going golden ratio:" << std::endl;

    for (auto eps : epsRange)
    {
        int m = 0;

        double L = b - a;
        double x1 = b - (b - a) / 1.618;
        double x2 = a + (b - a) / 1.618;

        double fx1 = f(x1);
        double fx2 = f(x2);

        do
        {
            if (fx1 < fx2)
            {
                b = x2;
                x2 = x1;
                fx2 = fx1;
                x1 = b - (b - a) / 1.618;
                fx1 = f(x1);
            }
            else
            {
                a = x1;
                x1 = x2;
                fx1 = fx2;
                x2 = a + (b - a) / 1.618;
                fx2 = f(x2);
            }

            L = b - a;
            m++;
        } while (L > eps);

        std::cout << "For eps = " << eps << " m = " << m << std::endl;
    }
}

int main()
{
    goldenRatio();
    return 0;
}*/


/*
#include <iostream>
#include <cmath>

// функция, определяющая уравнение, которое нужно решить
double f(double x) {
    return x + std::log(x);
}

// вычисление числа Фибоначчи по индексу
int fib(int n) {
    if (n <= 1) {
        return n;
    }
    else {
        return fib(n - 1) + fib(n - 2);
    }
}

void calcFibonacci() {
    double a = 0;
    double b = 10;
    double eps = 0.00001;

    int n = 0;
    while (fib(n) < (b - a) / eps) {
        n++;
    }

    int k = n - 1;
    double x1 = a + (double)fib(k - 1) / fib(k + 1) * (b - a);
    double x2 = a + (double)fib(k) / fib(k + 1) * (b - a);

    while (k >= 0) {
        if (f(x1) < f(x2)) {
            b = x2;
            x2 = x1;
            x1 = a + (double)fib(k - 1) / fib(k + 1) * (b - a);
        }
        else {
            a = x1;
            x1 = x2;
            x2 = a + (double)fib(k) / fib(k + 1) * (b - a);
        }

        k--;
    }

    double xmin = (a + b) / 2.0;
    double fmin = f(xmin);

    std::cout << "Minimum found: x = " << xmin << ", f(x) = " << fmin << std::endl;
}

int main() {
    calcFibonacci();
    return 0;
}*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//lab 01
/*
#include <iostream>
#include <cmath>
#include <chrono>

using namespace std;

double f(double x)
{
    return x + log(x);
}

// Метод дихотомии
double dichotomy(double a, double b, double eps)
{
    int m = 0;
    while (b - a > eps) {
        double x = (a + b) / 2.0;
        double fx = f(x);
        if (fx == 0.0) {
            return x;
        }
        if (f(a) * fx < 0.0) {
            b = x;
        }
        else {
            a = x;
        }
        m++;

        cout <<m<< " " << a<<" "<< b<< " "<< fx;
        cout << "Metod Dichotormei" << endl;    
        cout << "Peremenye:" << endl; 
        cout << "a = " << a << endl; 
        cout << "b = " << b << endl; 
        cout << "fx = " << fx << endl; 
        cout << "x = " << x << endl;
    }
    return (a + b) / 2.0;
    
}
     
// Метод золотого сечения
double goldenRatio(double a, double b, double eps)
{
    int m = 0;
    double phi = (1 + sqrt(5)) / 2.0;
    double x1 = b - (b - a) / phi;
    double x2 = a + (b - a) / phi;
    double fx1 = f(x1);
    double fx2 = f(x2);
    while (b - a > eps) {
        if (fx1 < fx2) {
            b = x2;
            x2 = x1;
            fx2 = fx1;
            x1 = b - (b - a) / phi;
            fx1 = f(x1);
        }
        else {
            a = x1;
            x1 = x2;
            fx1 = fx2;
            x2 = a + (b - a) / phi;
            fx2 = f(x2);
        }
        m++;
        
        cout << "" << endl;
        cout << "Metod Zolotoe sechenie" << endl;
        cout << "Peremenye:" << endl;
        cout << "a = " << a << endl;
        cout << "b = " << b << endl;
        cout << "fx1 = " << fx1 << endl;
        cout << "fx2 = " << fx2 << endl;
        cout << "x1 = " << x1 << endl;
        cout << "x2 = " << x2 << endl;

    }
    return (a + b) / 2.0;
}

// Метод Фибоначчи
double fibonacci(double a, double b, double eps)
{
    int m = 0;
    double L = b - a;
    double Fn1 = 1.0, Fn2 = 1.0;
    while (L / eps >= Fn2) {
        double x1 = a + (Fn1 / Fn2) * (b - a);
        double x2 = a + b - x1;
        double fx1 = f(x1);
        double fx2 = f(x2);
        if (fx1 <= fx2) {
            b = x2;
        }
        else {
            a = x1;
        }
        L = b - a;
        double temp = Fn2;
        Fn2 = Fn1 + Fn2;
        Fn1 = temp;
        m++;

        cout << "" << endl;
        cout << "Metod Fibonachi" << endl;
        cout << "Peremenye:" << endl;
        cout << "a = " << a << endl;
        cout << "b = " << b << endl;
        cout << "fx1 = " << fx1 << endl;
        cout << "fx2 = " << fx2 << endl;
        cout << "x1 = " << x1 << endl;
        cout << "x2 = " << x2 << endl;
        cout << "Fn1 = " << Fn1 << endl;
        cout << "Fn2 = " << Fn2 << endl;

    }
    return (a + b) / 2.0;
}

int main()
{
    double a = 2.0;
    double b = 20.0;
    double eps = 0.0001;
   
    // Для метода дихотомии
    auto start_time_d = chrono::high_resolution_clock::now();
    double res_d = dichotomy(a, b, eps);
    auto end_time_d = chrono::high_resolution_clock::now();
    auto time_d = chrono::duration_cast<chrono::microseconds>(end_time_d - start_time_d).count();


    // Для метода золотого сечения
    auto start_time_gs = chrono::high_resolution_clock::now();
    double res_gs = goldenRatio(a, b, eps);
    auto end_time_gs = chrono::high_resolution_clock::now();
    auto time_gs = chrono::duration_cast<chrono::microseconds>(end_time_gs - start_time_gs).count();

    // Для метода Фибоначчи
    auto start_time_f = chrono::high_resolution_clock::now();
    double res_f = fibonacci(a, b, eps);
    auto end_time_f = chrono::high_resolution_clock::now();
    auto time_f = chrono::duration_cast<chrono::microseconds>(end_time_f - start_time_f).count();

    // Вывод результатов и времени работы каждого метода
    cout << "Sravnenie metodov " <<"\n"<<endl;
    cout << "Result Dichotomy: " << res_d << ", Time: " << time_d << " microseconds" <<"\n"<< endl;
    cout << "Result Golden Section: " << res_gs << ", Time: " << time_gs << " microseconds" <<"\n"<< endl;
    cout << "Result Fibonacci: " << res_f << ", Time: " << time_f << " microseconds" <<"\n"<< endl;

    // Сравнение времени работы методов
    if (time_d < time_gs && time_d < time_f) {
        cout << "Dichotomy method is the fastest" <<"\n"<< endl;
    }
    else if (time_gs < time_d && time_gs < time_f) {
        cout << "Golden Section method is the fastest" <<"\n"<< endl;
    }
    else {
        cout << "Fibonacci method is the fastest" <<"\n"<< endl;
    }
}
*/
















//LAB1 3M
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Dihotormia 100%
/*
#include <iostream>
#include <cmath>
#include <chrono>

using namespace std;

double f(double x, double y) {
    return x * x + 2 * x + 2-y;
}

void dichotomy(double a, double b, double y, double eps) {
    int m = 0;
    double x1, x2, fx1, fx2;
    while (abs(b - a) > eps) {
        m++;
        x1 = (a + b) / 2.0;
        x2 = x1 + eps / 2.0;
        fx1 = f(x1, y);
        fx2 = f(x2, y);
        if (fx1 * fx2 > 0) {
            a = x1;
        }
        else {
            b = x2;
        }
        cout << "Iteration: " << m << ", a: " << a << ", b: " << b << ", x1: " << x1 << ", x2: " << x2 << ", f(x1): " << fx1 << ", f(x2): " << fx2 << endl;
    }
    double x = (a + b) / 2.0;
    double fx = f(x, y);
    cout << "Final result: x = " << x << ", f(x) = " << fx << endl;
}

int main() {
    double a = 2;
    double b = 20;
    double y = 5;
    double eps = 0.0001;

    auto start = std::chrono::high_resolution_clock::now();
    dichotomy(a, b, y, eps);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;

    return 0;
}

*/
//золотое сечения 100%
/*
#include <iostream>
#include <cmath>
#include <chrono>

using namespace std;

double f(double x, double y) {
    return x * x + 2 * x + 2 - y;
}

void golden_section(double a, double b, double y, double eps) {
    int m = 0;
    double phi = (1 + sqrt(5)) / 2;
    double resphi = 2 - phi;
    double x1 = b - resphi * (b - a);
    double x2 = a + resphi * (b - a);
    double fx1 = f(x1, y);
    double fx2 = f(x2, y);
    while ((b - a) / 2.0 > eps) {
        m++;
        cout << "Iteration: " << m << ", a: " << a << ", b: " << b << ", x1: " << x1 << ", x2: " << x2 << ", f(x1): " << fx1 << ", f(x2): " << fx2 << endl;
        if (fx1 < fx2) {
            b = x2;
            x2 = x1;
            fx2 = fx1;
            x1 = b - resphi * (b - a);
            fx1 = f(x1, y);
        }
        else {
            a = x1;
            x1 = x2;
            fx1 = fx2;
            x2 = a + resphi * (b - a);
            fx2 = f(x2, y);
        }
    }
    double x = (a + b) / 2.0;
    double fx = f(x, y);
    cout << "Iteration: " << m << ", a: " << a << ", b: " << b << ", x1: " << x1 << ", x2: " << x2 << ", f(x1): " << fx1 << ", f(x2): " << fx2 << endl;
   // cout << "  " << m << "  " << a << "  " << b << "  " << x1 << "  " << x2 << "  " << fx1 << ", " << fx2 << endl;

    cout << "Final result: x = " << x << ", f(x) = " << fx << endl;
}

int main() {
    double a = 2;
    double b = 20;
    double y = 5;
    double eps = 0.0001;
    golden_section(a, b, y, eps);
    auto start = std::chrono::high_resolution_clock::now();

    golden_section(a, b, y, eps);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;

    return 0;
}*/

//фибоначи 100%























/*
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>

using namespace std;

double f(double x, double y) {
    return x * x + 2 * x + 2 - y;
}

double fib(int n) {
    double phi = (1 + sqrt(5)) / 2;
    return (pow(phi, n) - pow(1 - phi, n)) / sqrt(5);
}

void fibonacci(double a, double b, double y, double eps) {
    int m = 0;
    int n = 0;
    while ((b - a) / fib(n + 2) > eps) {
        n++;
    }
    vector<double> d(n + 1);
    d[0] = (b - a) / fib(n + 2);
    d[1] = (b - a) / fib(n + 1);
    for (int i = 2; i <= n; i++) {
        d[i] = d[i - 1] / fib(i);
    }
    double x1 = a + d[n - 1];
    double x2 = b - d[n - 1];
    double fx1 = f(x1, y);
    double fx2 = f(x2, y);
    for (int i = 1; i < n - 1; i++) {
        if (fx1 > fx2) {
            a = x1;
            x1 = x2;
            x2 = b - d[n - i - 1];
            fx1 = fx2;
            fx2 = f(x2, y);
        }
        else {
            b = x2;
            x2 = x1;
            x1 = a + d[n - i - 1];
            fx2 = fx1;
            fx1 = f(x1, y);
        }
        m++;
        cout << "Iteration: " << m << ", a: " << a << ", b: " << b << ", x1: " << x1 << ", x2: " << x2 << ", f(x1): " << fx1 << ", f(x2): " << fx2 << endl;
       // cout << "  " << m << "  " << a << "  " << b << "  " << x1 << "  " << x2 << "  " << fx1 << ", " << fx2 << endl;
    }
    if (fx1 > fx2) {
        a = x1;
    }
    else {
        b = x2;
    }
    double x = (a + b) / 2.0;
    double fx = f(x, y);
    m++;
  //  cout << "Iteration: " << m << ", a: " << a << ", b: " << b << ", x1: " << x1 << ", x2: " << x2 << ", f(x1): " << fx1 << ", f(x2): " << fx2 << endl;
    cout << "  " << m << "  " << a << "  " << b << "  " << x1 << "  " << x2 << "  " << fx1 << ", " << fx2 << endl;
    cout << "Final result: x = " << x << ", f(x) = " << fx << endl;
}

int main() {
    double a = 2;
    double b = 20;
    double eps = 0.0001;
    double y = 5;
    fibonacci(a, b, y, eps);
    auto start = std::chrono::high_resolution_clock::now();

    fibonacci(a, b, y, eps);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;

    return 0;
}*/











///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
#include <iostream>
#include <cmath>
#include <chrono>

using namespace std;

double f(double x, double y) {
    return x + log(x) - y;
}

void dichotomy(double a, double b, double y, double eps) {
    int m = 0;
    double x1, x2, fx1, fx2;
    while (abs(b - a) > eps) {
        m++;
        x1 = (a + b) / 2.0;
        x2 = x1 + eps / 2.0;
        fx1 = f(x1, y);
        fx2 = f(x2, y);
        if (fx1 * fx2 > 0) {
            a = x1;
        }
        else {
            b = x2;
        }
        cout << "Iteration: " << m << ", a: " << a << ", b: " << b << ", x1: " << x1 << ", x2: " << x2 << ", f(x1): " << fx1 << ", f(x2): " << fx2 << endl;
    }
    double x = (a + b) / 2.0;
    double fx = f(x, y);
    cout << "Final result: x = " << x << ", f(x) = " << fx << endl;
}

int main() {
    double a = 2;
    double b = 20;
    double y = 5;
    double eps = 0.0001;

    auto start = std::chrono::high_resolution_clock::now();
    dichotomy(a, b, y, eps);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;

    return 0;
}*/
/*
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double func(double x, double y) {
    return x + log(x) - y;
}

vector<double> bisection(double a, double b, double y, double eps) {
    vector<double> results;
    double c = (a + b) / 2;
    while (fabs(b - a) > eps) {
        if (func(c, y) == 0.0) {
            results.push_back(c);
            return results;
        }
        else if (func(c, y) * func(a, y) < 0)
            b = c;
        else
            a = c;
        c = (a + b) / 2;
        results.push_back(c);
    }
    return results;
}

int main() {
    double y, eps;
    cout << "Enter y: ";
    cin >> y;
    cout << "Enter eps: ";
    cin >> eps;
    double a = exp(y) / (y + 1);
    double b = a * 2;
    vector<double> results = bisection(a, b, y, eps);
    cout << "Results:\n";
    cout << "x\tf(x)\n";
    for (double x : results)
        cout << x << "\t" << x + log(x) << "\n";
    return 0;
}

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double func(double x, double y) {
    return x + log(x) - y;
}

vector<double> dichotomy(double a, double b, double y, double eps) {
    vector<double> results;
    while (fabs(b - a) > eps) {
        double c = (a + b) / 2;
        if (func(c, y) == 0.0) {
            results.push_back(c);
            return results;
        }
        if (func(a, y) * func(c, y) < 0)
            b = c;
        else
            a = c;
        results.push_back(c);
    }
    return results;
}

int main() {
    double y, eps;
    cout << "Enter y: ";
    cin >> y;
    cout << "Enter eps: ";
    cin >> eps;
    double a = exp(y) / (y + 1);
    double b = a * 2;
    vector<double> results = dichotomy(a, b, y, eps);
    cout << "Results:\n";
    cout << "x\tf(x)\n";
    for (double x : results)
        cout << x << "\t" << x + log(x) << "\n";
    return 0;
}*/



//Золотое сечения минимум
/*
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double func(double x, double y) {
    return x + log(x) - y;
}

double goldenSection(double a, double b, double y, double eps) {
    double phi = (1 + sqrt(5)) / 2;
    double x1 = b - (b - a) / phi;
    double x2 = a + (b - a) / phi;
    double f1 = func(x1, y);
    double f2 = func(x2, y);
    vector<double> results;
    results.push_back(x1);
    results.push_back(x2);
    while (fabs(b - a) > eps) {
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - (b - a) / phi;
            f1 = func(x1, y);
        }
        else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + (b - a) / phi;
            f2 = func(x2, y);
        }
        results.push_back(x1);
    }
    cout << "Iterations:\n";
    cout << "x\tf(x)\n";
    for (double x : results)
        cout << x << "\t" << x + log(x) << "\n";
    return (a + b) / 2;
}

int main() {
    double y, eps;
    cout << "Enter y: ";
    cin >> y;
    cout << "Enter eps: ";
    cin >> eps;
    double a = exp(y) / (y + 1);
    double b = a * 2;
    double result = goldenSection(a, b, y, eps);
    cout << "Final result:\n";
    cout << "x\tf(x)\n";
    cout << result << "\t" << result + log(result) << "\n";
    return 0;
}*/

//Фибоначи минимум
/*
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double func(double x, double y) {
    return x + log(x) - y;
}

int fibonacci(int n) {
    if (n <= 0) return 0;
    if (n == 1) return 1;
    return fibonacci(n - 1) + fibonacci(n - 2);
}

double fibonacciSearch(double a, double b, double y, double eps) {
    int n = 0;
    while (fibonacci(n) < (b - a) / eps) n++;
    double x1 = a + (double)fibonacci(n - 2) / fibonacci(n) * (b - a);
    double x2 = a + (double)fibonacci(n - 1) / fibonacci(n) * (b - a);
    double f1 = func(x1, y);
    double f2 = func(x2, y);
    vector<double> results;
    results.push_back(x1);
    results.push_back(x2);
    for (int i = n - 3; i >= 0; i--) {
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (double)fibonacci(i) / fibonacci(i + 2) * (b - a);
            f1 = func(x1, y);
        }
        else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + (double)fibonacci(i + 1) / fibonacci(i + 2) * (b - a);
            f2 = func(x2, y);
        }
        results.push_back(x1);
    }
    cout << "Iterations:\n";
    cout << "x\tf(x)\n";
    for (int i = results.size() - 1; i >= 0; i--)
        cout << results[i] << "\t" << results[i] + log(results[i]) << "\n";
    return (a + b) / 2;
}

int main() {
    double y, eps;
    cout << "Enter y: ";
    cin >> y;
    cout << "Enter eps: ";
    cin >> eps;
    double a = exp(y) / (y + 1);
    double b = a * 2;
    double result = fibonacciSearch(a, b, y, eps);
    cout << "Final result:\n";
    cout << "x\tf(x)\n";
    cout << result << "\t" << result + log(result) << "\n";
    return 0;
}*/


//Дихотомия минимум
/*#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double func(double x, double y) {
    return x + log(x) - y;
}

void minimum(double a, double b, double y, double eps) {
    vector<pair<double, double>> results;
    while (fabs(b - a) > eps) {
        double c = (a + b) / 2;
        double fa = func(a, y);
        double fb = func(b, y);
        double fc = func(c, y);
        results.push_back(make_pair(c, fc));
        if (fa > fb) {
            a = c;
        }
        else {
            b = c;
        }
    }
    cout << "Iterations:\n";
    cout << "x\tf(x)\n";
    for (auto p : results)
        cout << p.first << "\t" << p.second << "\n";
}

int main() {
    double eps;
    cout << "Enter eps: ";
    cin >> eps;
    double a = 2;
    double b = 20;
    double y = 5;
    minimum(a, b, y, eps);
    return 0;
}*/

/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x, double y) {
    return x * x + 2 * x + 2 - y;
}

double golden_search(double y, double a, double b, double tol) {
    double phi = (1 + sqrt(5)) / 2;
    double x1 = b - (b - a) / phi;
    double x2 = a + (b - a) / phi;
    double fx1 = f(x1, y);
    double fx2 = f(x2, y);
    int iterations = 0;

    while (abs(b - a) > tol) {
        cout << "Iteration " << iterations << ": " << endl;
        cout << "x1 = " << x1 << ", x2 = " << x2 << ", a = " << a << ", b = " << b << ", f(x1) = " << fx1 << ", f(x2) = " << fx2 << endl;
        if (fx1 < fx2) {
            b = x2;
            x2 = x1;
            fx2 = fx1;
            x1 = b - (b - a) / phi;
            fx1 = f(x1, y);
        }
        else {
            a = x1;
            x1 = x2;
            fx1 = fx2;
            x2 = a + (b - a) / phi;
            fx2 = f(x2, y);
        }
        iterations++;
    }

    cout << "Minimum found at x = " << (a + b) / 2 << ", with a value of f(x) = " << f((a + b) / 2, y) << endl;
    return (a + b) / 2;
}

int main() {
    double y = 5;
    double a = 2;
    double b = 20;
    double tol = 0.0001;
    golden_search(y, a, b, tol);
    return 0;
}*/


/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return pow(x, 2) + 2 * x + 1;
}

double powell(double (*f)(double), double x0, double eps) {
    double x1 = x0, x2 = x0 + eps, x3 = x0 + 2 * eps;
    double fx1 = f(x1), fx2 = f(x2), fx3 = f(x3);
    double a, b;

    for (int i = 0; i < 100; i++) {
        double d1 = fx2 - fx1;
        double d2 = fx3 - fx2;
        double d3 = d2 - d1;

        a = (d3 * (x2 - x1) - d1 * (x3 - x2)) / (pow(x3 - x2, 2) * (x2 - x1) + pow(x2 - x1, 2) * (x3 - x2));
        b = (d2 - a * pow(x3 - x2, 2)) / (x3 - x2);

        double x4 = x3 + b * eps;
        double fx4 = f(x4);

        if (fx4 < fx3) {
            cout << "Iteration = " << i << ", x1 = " << x1 << ", x2 = " << x2 << ", x3 = " << x3 << ", fx1 = " << fx1 << ", fx2 = " << fx2 << ", fx3 = " << fx3 << ", a = " << a << ", b = " << b << endl;
            x1 = x2;
            x2 = x3;
            x3 = x4;
            fx1 = fx2;
            fx2 = fx3;
            fx3 = fx4;
        }
        else {
            eps = -eps / 2.0;
            x2 = x0 + eps;
            x3 = x0 + 2 * eps;
            fx2 = f(x2);
            fx3 = f(x3);
        }

        if (fabs(eps) < 1e-8) {
            break;
        }
    }

    return x2;
}

int main() {
    double x0 = 0.0;
    double eps = 0.01;

    double min_x = powell(f, x0, eps);

    cout << "Minimum point: " << min_x << endl;

    return 0;
}*/





/*
#include <iostream>
#include <cmath>

using namespace std;

// Define the function f(x) = x^2 + 2x + 1
double f(double x) {
    return pow(x, 2) + 2 * x + 1;
}

// Define the derivative of f(x), which is 2x + 2
double df(double x) {
    return 2 * x + 2;
}

// Implement the Newton method for finding the root of f(x)
double newton(double x0, double eps, int& iterations) {
    double x1 = x0 - f(x0) / df(x0);
    double fx1 = f(x1);
    double a = x0, b = x1;
    iterations = 1;
    cout << "x0 = " << x0 << ", f(x0) = " << f(x0) << endl;
    while (abs(fx1) > eps) {
        a = b;
        b = x1;
        x1 = x1 - f(x1) / df(x1);
        fx1 = f(x1);
        iterations++;
        cout << "x" << iterations << " = " << x1 << ", f(x" << iterations << ") = " << fx1 << endl;
    }
    cout << "Minimum point: " << x1 << endl;
    return x1;
}

// Implement the Newton-Raphson method for finding the root of f(x)
double newtonRaphson(double x0, double eps, int& iterations) {
    double x1 = x0 - f(x0) / df(x0);
    double fx1 = f(x1);
    double a = x0, b = x1;
    iterations = 1;
    cout << "x0 = " << x0 << ", f(x0) = " << f(x0) << endl;
    while (abs(x1 - a) > eps && abs(fx1) > eps) {
        a = b;
        b = x1;
        x1 = x1 - f(x1) / df(b);
        fx1 = f(x1);
        iterations++;
        cout << "x" << iterations << " = " << x1 << ", f(x" << iterations << ") = " << fx1 << endl;
    }
    cout << "Minimum point: " << x1 << endl;
    return x1;
}

int main() {
    double x0 = 0; // Initial guess
    double eps = 1e-6; // Tolerance for stopping criteria
    int iterations;

    cout << "Newton method:" << endl;
    double x_newton = newton(x0, eps, iterations);
    cout << "Number of iterations: " << iterations << endl;

    cout << endl;

    cout << "Newton-Raphson method:" << endl;
    double x_newtonRaphson = newtonRaphson(x0, eps, iterations);
    cout << "Number of iterations: " << iterations << endl;

    return 0;
}*/



/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return x * x + 2 * x + 1;
}

double df(double x) {
    return 2 * x + 2;
}

void middle_point_method(double a, double b, double eps) {
    double L = a;
    double R = b;
    double z;
    int i = 0;
    while (true) {
        z = (L + R) / 2;
        double df_z = df(z);
        cout << "Iteration " << i << ": z = " << z << ", f'(z) = " << df_z << endl;
        if (abs(df_z) <= eps) {
            break;
        }
        if (df_z < 0) {
            L = z;
        }
        else {
            R = z;
        }
        i++;
    }
    cout << "The minimum point is " << z << " with f(z) = " << f(z) << endl;
}

int main() {
    double a = -5;
    double b = 5;
    double eps = 0.01;
    middle_point_method(a, b, eps);
    return 0;
}*/



































//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//f(x) = x^2+2x+1


// ПАУЭЛА
/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return pow(x, 2) + 2 * x + 1;
}

double powell(double (*f)(double), double x0, double eps) {
    double x1 = x0, x2 = x0 + eps, x3 = x0 + 2 * eps;
    double fx1 = f(x1), fx2 = f(x2), fx3 = f(x3);
    double a, b;

    for (int i = 0; i < 100; i++) {
        double d1 = fx2 - fx1;
        double d2 = fx3 - fx2;
        double d3 = d2 - d1;

        a = (d3 * (x2 - x1) - d1 * (x3 - x2)) / (pow(x3 - x2, 2) * (x2 - x1) + pow(x2 - x1, 2) * (x3 - x2));
        b = (d2 - a * pow(x3 - x2, 2)) / (x3 - x2);

        double x4 = x3 + b * eps;
        double fx4 = f(x4);

        if (fx4 < fx3) {
            cout << "Iteration = " << i << ", x1 = " << x1 << ", x2 = " << x2 << ", x3 = " << x3 << ", fx1 = " << fx1 << ", fx2 = " << fx2 << ", fx3 = " << fx3 << ", a = " << a << ", b = " << b << endl;
            x1 = x2;
            x2 = x3;
            x3 = x4;
            fx1 = fx2;
            fx2 = fx3;
            fx3 = fx4;
        }
        else {
            eps = -eps / 2.0;
            x2 = x0 + eps;
            x3 = x0 + 2 * eps;
            fx2 = f(x2);
            fx3 = f(x3);
        }

        if (fabs(eps) < 1e-8) {
            break;
        }
    }

    return x2;
}

int main() {
    double x0 = 0.0;
    double eps = 0.01;

    double min_x = powell(f, x0, eps);

    cout << "Minimum point: " << min_x << endl;

    return 0;
}*/




//МЕТОД С ИСПОЛЬЗОВАНИЕМ ПРОИЗВОДНЫХ и МЕТОД НЬЮТОНА-РАФСОНА

/*
#include <iostream>
#include <cmath>

using namespace std;

// Определить функцию f(x) = x^2 + 2x + 1
double f(double x) {
    return pow(x, 2) + 2 * x + 1;
}

// Определите производную от f(x), которая равна 2x + 2
double df(double x) {
    return 2 * x + 2;
}

// Реализовать метод Ньютона для нахождения корня из f(x)
double newton(double x0, double eps, int& iterations) {
    double x1 = x0 - f(x0) / df(x0);
    double fx1 = f(x1);
    double a = x0, b = x1;
    iterations = 1;
    cout << "x0 = " << x0 << ", f(x0) = " << f(x0) << endl;
    while (abs(fx1) > eps) {
        a = b;
        b = x1;
        x1 = x1 - f(x1) / df(x1);
        fx1 = f(x1);
        iterations++;
        cout << "x" << iterations << " = " << x1 << ", f(x" << iterations << ") = " << fx1 << endl;
    }
    cout << "Minimum point: " << x1 << endl;
    return x1;
}

// Реализовать метод Ньютона-Рафсона для нахождения корня из f(x)
double newtonRaphson(double x0, double eps, int& iterations) {
    double x1 = x0 - f(x0) / df(x0);
    double fx1 = f(x1);
    double a = x0, b = x1;
    iterations = 1;
    cout << "x0 = " << x0 << ", f(x0) = " << f(x0) << endl;
    while (abs(x1 - a) > eps && abs(fx1) > eps) {
        a = b;
        b = x1;
        x1 = x1 - f(x1) / df(b);
        fx1 = f(x1);
        iterations++;
        cout << "x" << iterations << " = " << x1 << ", f(x" << iterations << ") = " << fx1 << endl;
    }
    cout << "Minimum point: " << x1 << endl;
    return x1;
}

int main() {
    double x0 = 0; // Первоначальное предположение
    double eps = 1e-6; 
    int iterations;

    cout << "Newton method:" << endl;
    double x_newton = newton(x0, eps, iterations);
    cout << "Number of iterations: " << iterations << endl;

    cout << endl;

    cout << "Newton-Raphson method:" << endl;
    double x_newtonRaphson = newtonRaphson(x0, eps, iterations);
    cout << "Number of iterations: " << iterations << endl;

    return 0;
}
*/




//МЕТОД СРЕДНЕЙ ТОЧКИ
/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return x * x + 2 * x + 1;
}

double df(double x) {
    return 2 * x + 2;
}

void middle_point_method(double a, double b, double eps) {
    double L = a;
    double R = b;
    double z;
    int i = 0;
    while (true) {
        z = (L + R) / 2;
        double df_z = df(z);
        cout << "Iteration " << i << ": z = " << z << ", f'(z) = " << df_z << endl;
        if (abs(df_z) <= eps) {
            break;
        }
        if (df_z < 0) {
            L = z;
        }
        else {
            R = z;
        }
        i++;
    }
    cout << "The minimum point is " << z << " with f(z) = " << f(z) << endl;
}

int main() {
    double a = -5;
    double b = 5;
    double eps = 0.01;
    middle_point_method(a, b, eps);
    return 0;
}*/





//МЕТОД СЕКУЩИХ 
/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return pow(x, 2) + 2 * x + 1;
}

double secant(double x0, double x1, double eps, int& iterations) {
    double fx0 = f(x0);
    double fx1 = f(x1);
    double x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
    double fx2 = f(x2);
    iterations = 1;
    cout << "x0 = " << x0 << ", f(x0) = " << fx0 << endl;
    cout << "x1 = " << x1 << ", f(x1) = " << fx1 << endl;
    while (abs(fx2) > eps) {
        x0 = x1;
        fx0 = fx1;
        x1 = x2;
        fx1 = fx2;
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
        fx2 = f(x2);
        iterations++;
        cout << "x" << iterations << " = " << x2 << ", f(x" << iterations << ") = " << fx2 << endl;
    }
    cout << "Minimum point: " << x2 << endl;
    return x2;
}

int main() {
    double x0 = 0; 
    double x1 = 1; 
    double eps = 1e-6; 
    int iterations;
    cout << "Secant method:" << endl;
    double x_secant = secant(x0, x1, eps, iterations);
    cout << "Number of iterations: " << iterations << endl;

    return 0;
}*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




















////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*15 Удалить из массива все элементы, расположенные между первым и
вторым нулевыми элементами. Если удаление элементов невозможно, выдать об этом сообщение.*/ 
/*
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main()
{
    const int size = 10;
    char arr[size];

    
    srand(time(NULL));
    for (int i = 0; i < size; i++)
    {
        arr[i] = rand() % 10 + '0'; 
        cout << arr[i] << " ";
    }
    cout << endl;

    int first_zero = -1, second_zero = -1;
    for (int i = 0; i < size; i++)
    {
        if (arr[i] == '0')
        {
            if (first_zero == -1)
            {
                first_zero = i;
            }
            else if (second_zero == -1)
            {
                second_zero = i;
                break;
            }
        }
    }

   
    if (first_zero == -1 || second_zero == -1)
    {
        cout << "Error: the array does not have two zero elements." << endl;
    }
    else
    {
        
        for (int i = 0; i <= first_zero; i++)
        {
            cout << arr[i] << " ";
        }
        for (int i = second_zero; i < size; i++)
        {
            cout << arr[i] << " ";
        }
        cout << endl;
    }

    return 0;
}
*/


/*15 После каждого элемента с отрицательным значением вставить элемент, равный абсолютной величине отрицательного элемента.
Если вставка элементов невозможна, выдать об этом сообщение.*/
/*
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

int main()
{
    const int size = 10;
    int arr[size * 2]; // увеличиваем размер массива в два раза

    // Заполнение массива случайными числами
    srand(time(NULL));
    for (int i = 0; i < size; i++)
    {
        arr[i] = rand() % 21 - 10; // генерируем случайное число в диапазоне [-10, 10]
        cout << arr[i] << " ";
    }
    cout << endl;

    // Вставка элементов после отрицательных элементов
    int new_size = size;
    for (int i = 0; i < new_size; i++)
    {
        if (arr[i] < 0)
        {
            // Определяем абсолютную величину отрицательного числа
            int abs_value = abs(arr[i]);
            // Увеличиваем размер массива на 1
            new_size++;
            // Сдвигаем элементы массива справа от текущего на одну позицию вправо
            for (int j = new_size - 1; j > i + 1; j--)
            {
                arr[j] = arr[j - 1];
            }
            // Вставляем новый элемент после текущего
            arr[i + 1] = abs_value;
            // Увеличиваем индекс текущего элемента на 1, чтобы не обрабатывать вставленный элемент
            i++;
        }
    }

    // Проверка наличия вставленных элементов
    if (new_size == size)
    {
        cout << "Error: the array does not have negative elements." << endl;
    }
    else if (new_size > size * 2) // проверяем, что новый размер массива не превышает его выделенный размер
    {
        cout << "Error: cannot insert more elements." << endl;
    }
    else
    {
        // Вывод результата
        for (int i = 0; i < new_size; i++)
        {
            cout << arr[i] << " ";
        }
        cout << endl;
    }

    return 0;
}*/




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*СОРТИРОВКА[3, 7, 24] Сортировка вставками(Insertion Sort) Сортировка подсчётом(Counting Sort)
Сортировка змейкой (Snake Sort)*/

/*
#include <iostream>
#include <vector>

void insertionSort(std::vector<int>& arr) {
    int n = arr.size();

    for (int i = 1; i < n; i++) {
        int key = arr[i];
        int j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }

        arr[j + 1] = key;

        for (int k = 0; k < n; k++) {
            std::cout << arr[k] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    std::vector<int> arr = { 9, 1, 8, 2, 7, 3, 6, 4, 5, 10 };
    std::cout << "Unsorted array: ";
    for (int i = 0; i < arr.size(); i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
    insertionSort(arr);
    std::cout << "Sorted array: ";
    for (int i = 0; i < arr.size(); i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
    return 0;
}*/

/*
#include <iostream>
#include <vector>

void countingSort(std::vector<int>& arr)
{
    // Находим максимальное и минимальное значения в массиве
    int maxVal = arr[0], minVal = arr[0];
    for (int i = 1; i < arr.size(); ++i) {
        if (arr[i] > maxVal) maxVal = arr[i];
        if (arr[i] < minVal) minVal = arr[i];
    }

    // Создаем массив для подсчета количества элементов
    std::vector<int> count(maxVal - minVal + 1, 0);

    // Считаем количество элементов каждого значения в массиве
    for (int i = 0; i < arr.size(); ++i) {
        count[arr[i] - minVal]++;
    }

    // Пересчитываем элементы в массиве с учетом количества
    int pos = 0;
    for (int i = 0; i < count.size(); ++i) {
        for (int j = 0; j < count[i]; ++j) {
            arr[pos++] = i + minVal;

            // Выводим текущее состояние массива
            std::cout << "Iteration " << pos << ": ";
            for (int k = 0; k < arr.size(); ++k) {
                std::cout << arr[k] << " ";
            }
            std::cout << std::endl;
        }
    }
}

int main()
{
    std::vector<int> arr = { 9, 1, 8, 2, 7, 3, 6, 4, 5, 10 };

    // Выводим исходное состояние массива
    std::cout << "Iteration 0: ";
    for (int i = 0; i < arr.size(); ++i) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;

    // Сортируем массив
    countingSort(arr);

    return 0;
}*/


/*
#include <iostream>
#include <vector>
#include <algorithm>

void snakeSort(std::vector<int>& arr) {
    int n = arr.size();
    bool swapped = true;
    int start = 0, end = n - 1;

    while (swapped) {
        swapped = false;

        for (int i = start; i < end; i++) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }

        if (!swapped) {
            break;
        }

        swapped = false;

        for (int i = end - 1; i >= start; i--) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }

        start++;
        end--;

        for (int i = 0; i < n; i++) {
            std::cout << arr[i] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    std::vector<int> arr = { 9, 1, 8, 2, 7, 3, 6, 4, 5, 10 };
    std::cout << "Unsorted array: ";
    for (int i = 0; i < arr.size(); i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
    snakeSort(arr);
    std::cout << "Sorted array: ";
    for (int i = 0; i < arr.size(); i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
    return 0;
}*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////ОБЩИЙ ДЛЯ 3-Х
/*
#include <iostream>
#include <vector>
#include <algorithm>

void insertionSort(std::vector<int>& arr) {
    int n = arr.size();

    for (int i = 1; i < n; i++) {
        int key = arr[i];
        int j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }

        arr[j + 1] = key;

        std::cout << "Insertion sort iteration " << i << ": ";
        for (int k = 0; k < n; k++) {
            std::cout << arr[k] << " ";
        }
        std::cout << std::endl;
    }
}

void countingSort(std::vector<int>& arr)
{
    int maxVal = arr[0], minVal = arr[0];
    for (int i = 1; i < arr.size(); ++i) {
        if (arr[i] > maxVal) maxVal = arr[i];
        if (arr[i] < minVal) minVal = arr[i];
    }

    std::vector<int> count(maxVal - minVal + 1, 0);

    for (int i = 0; i < arr.size(); ++i) {
        count[arr[i] - minVal]++;
    }

    int pos = 0;
    for (int i = 0; i < count.size(); ++i) {
        for (int j = 0; j < count[i]; ++j) {
            arr[pos++] = i + minVal;

            std::cout << "Counting sort iteration " << pos << ": ";
            for (int k = 0; k < arr.size(); ++k) {
                std::cout << arr[k] << " ";
            }
            std::cout << std::endl;
        }
    }
}

void snakeSort(std::vector<int>& arr) {
    int n = arr.size();
    bool swapped = true;
    int start = 0, end = n - 1;

    while (swapped) {
        swapped = false;

        for (int i = start; i < end; i++) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }

        if (!swapped) {
            break;
        }

        swapped = false;

        for (int i = end - 1; i >= start; i--) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }

        start++;
        end--;

        std::cout << "Snake sort iteration " << start << ": ";
        for (int i = 0; i < n; i++) {
            std::cout << arr[i] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    std::vector<int> arr = { 9, 1, 8, 2, 7, 3, 6, 4, 5, 10 };

    std::cout << "Unsorted array: ";
    for (int i = 0; i < arr.size(); i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;



std::vector<int> arr2 = arr;
std::cout << "Sorting using insertion sort:\n";
insertionSort(arr2);

std::vector<int> arr3 = arr;
std::cout << "\nSorting using counting sort:\n";
countingSort(arr3);

std::vector<int> arr4 = arr;
std::cout << "\nSorting using snake sort:\n";
snakeSort(arr4);

return 0;
}*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





/*
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main()
{
    const int size = 10;
    char arr[size];

    // Заполнение массива случайными числами
    srand(time(NULL));
    for (int i = 0; i < size; i++)
    {
        arr[i] = rand() % 10 + '0'; // Преобразуем случайное число в символ
        cout << arr[i] << " ";
    }
    cout << endl;

    // Поиск первого и второго нулевых элементов
    int first_zero = -1, second_zero = -1;
    for (int i = 0; i < size; i++)
    {
        if (arr[i] == '0')
        {
            if (first_zero == -1)
            {
                first_zero = i;
            }
            else if (second_zero == -1)
            {
                second_zero = i;
                break;
            }
        }
    }

    // Проверка наличия двух нулевых элементов
    if (first_zero == -1 || second_zero == -1)
    {
        cout << "Error: the array does not have two zero elements." << endl;
    }
    else
    {
        // Вывод измененного массива без элементов между первым и вторым нулевыми элементами
        for (int i = 0; i <= first_zero; i++)
        {
            cout << arr[i] << " ";
        }
        for (int i = second_zero; i < size; i++)
        {
            cout << arr[i] << " ";
        }
        cout << endl;

        // Проверка наличия нулей между первым и вторым нулевыми элементами
        if (second_zero - first_zero > 1)
        {
            // Вывод нулей между первым и вторым нулевыми элементами
            for (int i = first_zero + 1; i < second_zero; i++)
            {
                cout << '0' << " ";
            }
            cout << endl;
        }
    }

    return 0;
}*/













/*
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main()
{
    const int size = 10;
    char arr[size];

    // Заполнение массива случайными числами
    srand(time(NULL));
    for (int i = 0; i < size; i++)
    {
        arr[i] = rand() % 10 + '0'; // Преобразуем случайное число в символ
        cout << arr[i] << " ";
    }
    cout << endl;

    // Поиск первого и второго нулевых элементов
    int first_zero = -1, second_zero = -1;
    for (int i = 0; i < size; i++)
    {
        if (arr[i] == '0')
        {
            if (first_zero == -1)
            {
                first_zero = i;
            }
            else if (second_zero == -1)
            {
                second_zero = i;
                break;
            }
        }
    }

    // Проверка наличия двух нулевых элементов
    if (first_zero == -1 || second_zero == -1)
    {
        cout << "Error: the array does not have two zero elements." << endl;
    }
    else
    {
        // Вывод измененного массива без элементов между первым и вторым нулевыми элементами
        for (int i = 0; i <= first_zero; i++)
        {
            cout << arr[i] << " ";
        }
        for (int i = second_zero; i < size; i++)
        {
            cout << arr[i] << " ";
        }
        cout << endl;
    }

    return 0;
}*/










/*
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;

int main()
{
    const int size = 10;
    int arr[size * 2]; // увеличиваем размер массива в два раза

    // Заполнение массива случайными числами
    srand(time(NULL));
    for (int i = 0; i < size; i++)
    {
        arr[i] = rand() % 21 - 10; // генерируем случайное число в диапазоне [-10, 10]
        cout << arr[i] << " ";
    }
    cout << endl;

    // Вставка элементов после отрицательных элементов
    int new_size = size;
    for (int i = 0; i < new_size; i++)
    {
        if (arr[i] < 0)
        {
            // Определяем абсолютную величину отрицательного числа
            int abs_value = abs(arr[i]);
            // Увеличиваем размер массива на 1
            new_size++;
            // Сдвигаем элементы массива справа от текущего на одну позицию вправо
            for (int j = new_size - 1; j > i + 1; j--)
            {
                arr[j] = arr[j - 1];
            }
            // Вставляем новый элемент после текущего
            arr[i + 1] = abs_value;
            // Увеличиваем индекс текущего элемента на 1, чтобы не обрабатывать вставленный элемент
            i++;
        }
    }

    // Проверка наличия вставленных элементов
    if (new_size == size)
    {
        cout << "Error: the array does not have negative elements." << endl;
    }
    else if (new_size > size * 2) // проверяем, что новый размер массива не превышает его выделенный размер
    {
        cout << "Error: cannot insert more elements." << endl;
    }
    else
    {
        // Вывод результата
        for (int i = 0; i < new_size; i++)
        {
            cout << arr[i] << " ";
        }
        cout << endl;
    }

    return 0;
}*/
































///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//1 Задание
/*
#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream input("C:/Users/ASUS/Desktop/books.txt");
    if (!input.is_open()) {
        std::cerr << "Failed to open file\n";
        return 1;
    }

    std::string line;
    int count = 0;
    while (std::getline(input, line) && count < 3) {
        std::cout << line << std::endl;
        count++;
    }

    input.close();
    return 0;
}*/


//2 Задание 
/*
#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream input("C:/Users/ASUS/Desktop/series.txt");
    if (!input.is_open()) {
        std::cerr << "Failed to open file\n";
        return 1;
    }

    std::string line;
    int num_lines = 0;
    int num_words = 0;
    int num_chars = 0;
    while (std::getline(input, line)) {
        num_lines++;

       
        int line_words = 0;
        bool in_word = false;
        for (char c : line) {
            if (isspace(c)) {
                if (in_word) {
                    line_words++;
                    in_word = false;
                }
            }
            else {
                in_word = true;
            }
        }
        if (in_word) {
            line_words++;
        }
        num_words += line_words;

       
        for (char c : line) {
            if (c != ' ' && c != '.') {
                num_chars++;
            }
        }
    }

    std::cout << "kol-vo strok: " << num_lines << std::endl;
    std::cout << "kol-vo slov: " << num_words << std::endl;
    std::cout << "kol-vo symvolov: " << num_chars << std::endl;

    input.close();
    return 0;
}*/


//3 Задание
/*
#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream input("C:/Users/ASUS/Desktop/books.txt");
    if (!input.is_open()) {
        std::cerr << "Failed to open file\n";
        return 1;
    }

    std::string line;
    std::string longest_word;
    while (std::getline(input, line)) {
        std::string word;
        for (char c : line) {
            if (c != ' ' && c != ',' && c != '.' && c != '-' && c != '(' && c != ')') {
                word += c;
            }
            else {
                if (word.length() > longest_word.length()) {
                    longest_word = word;
                }
                word = "";
            }
        }
        if (word.length() > longest_word.length()) {
            longest_word = word;
        }
    }

    std::cout << "Samoe dlinnoe slovo: " << longest_word << std::endl;

    input.close();
    return 0;
}*/


//4 Задание 
/*
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

using namespace std;

int main()
{
    unordered_map<string, int> fruit_counts;
    string fruit_name;
    ifstream fruit_file("C:/Users/ASUS/Desktop/fruit.txt");

    while (fruit_file >> fruit_name) {
        ++fruit_counts[fruit_name];
    }

    cout << "Nazvanya vstrechaytsa:\n";
    for (const auto& pair : fruit_counts) {
        cout << "\"" << pair.first << "\" - " << pair.second << " raz\n";
    }

    return 0;
}*/


//5 Задание 
/*
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
int main() {
    std::ifstream firstFile("C:/Users/ASUS/Desktop/first.txt");
    std::ifstream secondFile("C:/Users/ASUS/Desktop/second.txt");

    if (!firstFile.is_open() || !secondFile.is_open()) {
        std::cout << "Не удалось открыть файлы!" << std::endl;
        return 1;
    }

    std::vector<std::string> names;
    std::vector<std::string> positions;

    std::string name;
    while (std::getline(firstFile, name)) {
        names.push_back(name);
    }

    std::string position;
    while (std::getline(secondFile, position)) {
        positions.push_back(position);
    }

    if (names.size() != positions.size()) {
        std::cout << "Некорректные данные в файлах!" << std::endl;
        return 1;
    }

    for (size_t i = 0; i < names.size(); i++) {
        std::cout << "Sotrydnik: " << names[i] << ", dolshnost - " << positions[i] << std::endl;
    }

    firstFile.close();
    secondFile.close();

    return 0;
}
*/


//6 Задание
/*
#include <iostream>
#include <fstream>

int main() {
    std::ofstream outFile("alphabet.txt");

    if (!outFile.is_open()) {
        std::cout << "Не удалось открыть файл для записи!" << std::endl;
        return 1;
    }

    for (char c = 'A'; c <= 'Z'; c++) {
        outFile << c << std::endl;
        outFile << (char)(c + 32) << std::endl;
    }

    outFile.close();

    return 0;
}*/


//7 Задание 
/*
#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream incomeFile("C:/Users/ASUS/Desktop/income.txt");
    std::ifstream outcomeFile("C:/Users/ASUS/Desktop/outcome.txt");

    if (!incomeFile.is_open() || !outcomeFile.is_open()) {
        std::cout << "Не удалось открыть файлы!" << std::endl;
        return 1;
    }

    int incomeSum = 0;
    std::string incomeLine;
    while (std::getline(incomeFile, incomeLine)) {
        int income = std::stoi(incomeLine.substr(3)); 
        incomeSum += income;
    }

    int outcomeSum = 0;
    std::string outcomeLine;
    while (std::getline(outcomeFile, outcomeLine)) {
        int outcome = std::stoi(outcomeLine.substr(4)); 
        outcomeSum += outcome;
    }

    int profit = incomeSum + outcomeSum;
    std::cout << "Pribl: " << profit << " RUB" << std::endl;

    incomeFile.close();
    outcomeFile.close();

    return 0;
}
*/


//8 Задание 
/*
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>

int main() {
    std::ifstream file1("C:/Users/ASUS/Desktop/grades.txt");
    std::string line;
    std::vector<std::pair<std::string, double>> result;

    while (std::getline(file1, line)) {
        std::istringstream iss(line);
        std::vector<int> grades;
        std::string surname, name;
        iss >> surname >> name;
        int grade;
        while (iss >> grade) {
            grades.push_back(grade);
        }
        double avg_grade = static_cast<double>(std::accumulate(grades.begin(), grades.end(), 0)) / grades.size();
        if (avg_grade >= 4.5) {
            result.emplace_back(surname + " " + name.substr(0, 1) + ".", avg_grade);
        }
    }
    file1.close();

    for (const auto& pair : result) {
        std::cout << pair.first << ", sredny: " << std::fixed << std::setprecision(2) << pair.second << std::endl;
    }

    return 0;
}*/




//9 Задание 
/*
#include <iostream>
#include <fstream>
#include <map>
#include <string>

using namespace std;

int main() {
   
    map<char, string> dict = {
        {'q', "Q"}, {'w', "W"}, {'e', "E"}, {'r', "R"}, {'t', "T"}, {'y', "Y"}, {'u', "U"}, {'i', "I"}, {'o', "O"}, {'p', "P"},
        {'a', "A"}, {'s', "S"}, {'d', "D"}, {'f', "F"}, {'g', "G"}, {'h', "H"}, {'j', "J"}, {'k', "K"}, {'l', "L"}, {'z', "Z"},
        {'x', "X"}, {'c', "C"}, {'v', "V"}, {'b', "B"}, {'n', "N"}, {'m', "M"}, {',', ","}, {'.', "."}, {'?', "?"}, {'!', "!"},
    };

    
    string input;
    cout << "Stroka: ";
    getline(cin, input);

   
    string output;
    for (char c : input) {
        auto it = dict.find(c);
        if (it != dict.end()) {
            output += it->second;
        }
        else {
            output += c;
        }
    }

    
    ofstream out("C:/Users/ASUS/Desktop/result.txt");
    out << output;
    out.close();

    return 0;
}
*/


//10 Задание
/*
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


int timeDifference(string startTime, string endTime) {
    int startHour = stoi(startTime.substr(0, 2));
    int startMinute = stoi(startTime.substr(3, 2));
    int endHour = stoi(endTime.substr(0, 2));
    int endMinute = stoi(endTime.substr(3, 2));
    int totalTime = (endHour - startHour) * 60 + (endMinute - startMinute);
    return totalTime;
}

int main() {
    string logFileName = "C:/Users/ASUS/Desktop/crm_log.txt";
    string outputFileName = "C:/Users/ASUS/Desktop/best_employees.txt";
    ifstream logFile(logFileName);
    ofstream outputFile(outputFileName);
    string line;
    while (getline(logFile, line)) {
        
        string name = line.substr(0, line.find(","));
        string startTime = line.substr(line.find(",") + 2, 5);
        string endTime = line.substr(line.find(",") + 9, 5);
      
        int timeDiff = timeDifference(startTime, endTime);
     
        if (timeDiff > 240) {
            outputFile << name << endl;
        }
    }
 
    logFile.close();
    outputFile.close();
    return 0;
}*/
  




/*
#include <iostream>
#include <fstream>
#include <vector>

void insertionSort(std::vector<int>& arr) {
    int n = arr.size();

    for (int i = 1; i < n; i++) {
        int key = arr[i];
        int j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }

        arr[j + 1] = key;

        for (int k = 0; k < n; k++) {
            std::cout << arr[k] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    std::ifstream inFile("C:\\Users\\ASUS\\Desktop\\i.txt");
    if (!inFile) {
        std::cerr << "Unable to open file" << std::endl;
        return 1;
    }

    std::vector<int> arr;
    int num;
    while (inFile >> num) {
        arr.push_back(num);
    }
    inFile.close();

    std::cout << "Unsorted array: ";
    for (int i = 0; i < arr.size(); i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;

    insertionSort(arr);

    std::cout << "Sorted array: ";
    for (int i = 0; i < arr.size(); i++) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
*/




/*
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
void insertionSort(std::vector<int>& arr) {
    int n = arr.size();

    for (int i = 1; i < n; i++) {
        int key = arr[i];
        int j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }

        arr[j + 1] = key;

        std::cout << "Insertion sort iteration " << i << ": ";
        for (int k = 0; k < n; k++) {
            std::cout << arr[k] << " ";
        }
        std::cout << std::endl;
    }
}

void countingSort(std::vector<int>& arr)
{
    int maxVal = arr[0], minVal = arr[0];
    for (int i = 1; i < arr.size(); ++i) {
        if (arr[i] > maxVal) maxVal = arr[i];
        if (arr[i] < minVal) minVal = arr[i];
    }

    std::vector<int> count(maxVal - minVal + 1, 0);

    for (int i = 0; i < arr.size(); ++i) {
        count[arr[i] - minVal]++;
    }

    int pos = 0;
    for (int i = 0; i < count.size(); ++i) {
        for (int j = 0; j < count[i]; ++j) {
            arr[pos++] = i + minVal;

            std::cout << "Counting sort iteration " << pos << ": ";
            for (int k = 0; k < arr.size(); ++k) {
                std::cout << arr[k] << " ";
            }
            std::cout << std::endl;
        }
    }
}

void snakeSort(std::vector<int>& arr) {
    int n = arr.size();
    bool swapped = true;
    int start = 0, end = n - 1;

    while (swapped) {
        swapped = false;

        for (int i = start; i < end; i++) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }

        if (!swapped) {
            break;
        }

        swapped = false;

        for (int i = end - 1; i >= start; i--) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }

        start++;
        end--;

        std::cout << "Snake sort iteration " << start << ": ";
        for (int i = 0; i < n; i++) {
            std::cout << arr[i] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {

    std::ifstream infile("C:\\Users\\ASUS\\Desktop\\i.txt");
    std::vector<int> arr;

    int num;
    while (infile >> num) {
        arr.push_back(num);
    }
    std::cout << std::endl;



    std::vector<int> arr2 = arr;
    std::cout << "Sorting using insertion sort:\n";
    insertionSort(arr2);

    std::vector<int> arr3 = arr;
    std::cout << "\nSorting using counting sort:\n";
    countingSort(arr3);

    std::vector<int> arr4 = arr;
    std::cout << "\nSorting using snake sort:\n";
    snakeSort(arr4);

    return 0;
}*/








///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// из файла в файл 1
/*
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

void insertionSort(std::vector<int>& arr, std::ofstream& outfile) {
    int n = arr.size();

    for (int i = 1; i < n; i++) {
        int key = arr[i];
        int j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }

        arr[j + 1] = key;

        outfile << "Insertion sort iteration " << i << ": ";
        for (int k = 0; k < n; k++) {
            outfile << arr[k] << " ";
        }
        outfile << std::endl;
    }
}

void countingSort(std::vector<int>& arr, std::ofstream& outfile)
{
    int maxVal = arr[0], minVal = arr[0];
    for (int i = 1; i < arr.size(); ++i) {
        if (arr[i] > maxVal) maxVal = arr[i];
        if (arr[i] < minVal) minVal = arr[i];
    }

    std::vector<int> count(maxVal - minVal + 1, 0);

    for (int i = 0; i < arr.size(); ++i) {
        count[arr[i] - minVal]++;
    }

    int pos = 0;
    for (int i = 0; i < count.size(); ++i) {
        for (int j = 0; j < count[i]; ++j) {
            arr[pos++] = i + minVal;

            outfile << "Counting sort iteration " << pos << ": ";
            for (int k = 0; k < arr.size(); ++k) {
                outfile << arr[k] << " ";
            }
            outfile << std::endl;
        }
    }
}

void snakeSort(std::vector<int>& arr, std::ofstream& outfile) {
    int n = arr.size();
    bool swapped = true;
    int start = 0, end = n - 1;

    while (swapped) {
        swapped = false;

        for (int i = start; i < end; i++) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }

        if (!swapped) {
            break;
        }

        swapped = false;

        for (int i = end - 1; i >= start; i--) {
            if (arr[i] > arr[i + 1]) {
                std::swap(arr[i], arr[i + 1]);
                swapped = true;
            }
        }

        start++;
        end--;

        outfile << "Snake sort iteration " << start << ": ";
        for (int i = 0; i < n; i++) {
            outfile << arr[i] << " ";
        }
        outfile << std::endl;
    }
}

int main() {
    std::ifstream infile("C:\\Users\\ASUS\\Desktop\\i.txt");
    std::vector<int> arr;

    int num;
    while (infile >> num) {
        arr.push_back(num);
    }

    std::ofstream outfile("C:\\Users\\ASUS\\Desktop\\algor1.txt");

    outfile << "Input array: ";
    for (int i = 0; i < arr.size(); i++) {
        outfile << arr[i] << " ";
    }
    outfile << std::endl << std::endl;

    std::vector<int> arr1 = arr;
    insertionSort(arr1, outfile);
    outfile << std::endl;

    std::vector<int> arr2 = arr;
    countingSort(arr2, outfile);
    outfile << std::endl;

    std::vector<int> arr3 = arr;
    snakeSort(arr3, outfile);
    outfile << std::endl;

    outfile.close();

    return 0;
}*/




//2 Задание
/*
#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream input("C:\\Users\\ASUS\\Desktop\\algor.txt");
    if (!input.is_open()) {
        std::cerr << "Failed to open file\n";
        return 1;
    }

    std::string line;
    int count = 0;
    while (std::getline(input, line) && count < 3) {
        std::cout << line << std::endl;
        count++;
    }

    input.close();
    return 0;
}*/


//3 Задание 
/*
#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream input("C:\\Users\\ASUS\\Desktop\\algor.txt");
    if (!input.is_open()) {
        std::cerr << "Failed to open file\n";
        return 1;
    }

    std::string line;
    int num_lines = 0;
    int num_words = 0;
    int num_chars = 0;
    while (std::getline(input, line)) {
        num_lines++;


        int line_words = 0;
        bool in_word = false;
        for (char c : line) {
            if (isspace(c)) {
                if (in_word) {
                    line_words++;
                    in_word = false;
                }
            }
            else {
                in_word = true;
            }
        }
        if (in_word) {
            line_words++;
        }
        num_words += line_words;


        for (char c : line) {
            if (c != ' ' && c != '.') {
                num_chars++;
            }
        }
    }

    std::cout << "kol-vo strok: " << num_lines << std::endl;
    std::cout << "kol-vo slov: " << num_words << std::endl;
    std::cout << "kol-vo symvolov: " << num_chars << std::endl;

    input.close();
    return 0;
}
*/

//4 Задание
/*
#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream input("C:\\Users\\ASUS\\Desktop\\algor.txt");
    if (!input.is_open()) {
        std::cerr << "Failed to open file\n";
        return 1;
    }

    std::string line;
    std::string longest_word;
    while (std::getline(input, line)) {
        std::string word;
        for (char c : line) {
            if (c != ' ' && c != ',' && c != '.' && c != '-' && c != '(' && c != ')') {
                word += c;
            }
            else {
                if (word.length() > longest_word.length()) {
                    longest_word = word;
                }
                word = "";
            }
        }
        if (word.length() > longest_word.length()) {
            longest_word = word;
        }
    }

    std::cout << "Samoe dlinnoe slovo: " << longest_word << std::endl;

    input.close();
    return 0;
}*/


//5 Задание 
/*
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

using namespace std;

int main()
{
    unordered_map<string, int> fruit_counts;
    string fruit_name;
    ifstream fruit_file("C:\\Users\\ASUS\\Desktop\\algor.txt");

    while (fruit_file >> fruit_name) {
        ++fruit_counts[fruit_name];
    }

    cout << "Nazvanya vstrechaytsa:\n";
    for (const auto& pair : fruit_counts) {
        cout << "\"" << pair.first << "\" - " << pair.second << " raz\n";
    }

    return 0;
}
*/

//6 Задание 
/*
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main() {
    std::ifstream firstFile("C:\\Users\\ASUS\\Desktop\\algor.txt");
    std::ifstream secondFile("C:\\Users\\ASUS\\Desktop\\algorv.txt");

    if (!firstFile.is_open() || !secondFile.is_open()) {
        std::cout << "Не удалось открыть файлы!" << std::endl;
        return 1;
    }

    std::vector<std::string> fileContents;

    std::string line1;
    std::string line2;

    while (std::getline(firstFile, line1) && std::getline(secondFile, line2)) {
        fileContents.push_back(line1 + " " + line2);
    }

    if (fileContents.empty()) {
        std::cout << "Файлы пусты!" << std::endl;
        return 1;
    }

    for (const auto& line : fileContents) {
        std::cout << line << std::endl;
    }

    firstFile.close();
    secondFile.close();

    return 0;
}*/




//7 Задание
/*
#include <iostream>
#include <fstream>

int main() {
    std::ofstream outFile("C:\\Users\\ASUS\\Desktop\\n.txt");

    if (!outFile.is_open()) {
        std::cout << "Не удалось открыть файл для записи!" << std::endl;
        return 1;
    }

    for (int i = 1; i <= 36; i++) {
        outFile << i << std::endl;
    }

    outFile.close();

    return 0;
}*/


//8 Задание 
/*
#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream inFile("C:\\Users\\ASUS\\Desktop\\algor.txt");

    if (!inFile.is_open()) {
        std::cout << "Не удалось открыть файл для чтения!" << std::endl;
        return 1;
    }

    int insertionCount = 1, countingCount = 0, snakeCount = 0;
    std::string line;

    while (std::getline(inFile, line)) {
        if (line.find("Insertion sort iteration") != std::string::npos) {
            insertionCount++;
        }
        else if (line.find("Counting sort iteration") != std::string::npos) {
            countingCount++;
        }
        else if (line.find("Snake sort iteration") != std::string::npos) {
            snakeCount++;
        }
    }

    std::cout << "Insertions sort iteration = " << insertionCount << std::endl;
    std::cout << "Counting sort iteration = " << countingCount << std::endl;
    std::cout << "Snake sort iteration = " << snakeCount << std::endl;

    inFile.close();

    return 0;
}*/



//9 Задание 
/*
#include <iostream>
#include <fstream>
#include <map>
#include <string>

using namespace std;

int main() {

    map<char, string> dict = {
        {'q', "Q"}, {'w', "W"}, {'e', "E"}, {'r', "R"}, {'t', "T"}, {'y', "Y"}, {'u', "U"}, {'i', "I"}, {'o', "O"}, {'p', "P"},
        {'a', "A"}, {'s', "S"}, {'d', "D"}, {'f', "F"}, {'g', "G"}, {'h', "H"}, {'j', "J"}, {'k', "K"}, {'l', "L"}, {'z', "Z"},
        {'x', "X"}, {'c', "C"}, {'v', "V"}, {'b', "B"}, {'n', "N"}, {'m', "M"}, {',', ","}, {'.', "."}, {'?', "?"}, {'!', "!"},
    };


    string input;
    cout << "Stroka: ";
    getline(cin, input);


    string output;
    for (char c : input) {
        auto it = dict.find(c);
        if (it != dict.end()) {
            output += it->second;
        }
        else {
            output += c;
        }
    }


    ofstream out("C:/Users/ASUS/Desktop/result.txt");
    out << output;
    out.close();

    return 0;
}*/



//10 Задание
/*
#include <iostream>
#include <fstream>
#include <string>

int main() {
    std::ifstream inFile("C:\\Users\\ASUS\\Desktop\\algor.txt");

    if (!inFile.is_open()) {
        std::cout << "Не удалось открыть файл для чтения!" << std::endl;
        return 1;
    }

    int insertionCount = 1, countingCount = 0, snakeCount = 0;
    std::string line;

    while (std::getline(inFile, line)) {
        if (line.find("Insertion sort iteration") != std::string::npos) {
            insertionCount++;
        }
        else if (line.find("Counting sort iteration") != std::string::npos) {
            countingCount++;
        }
        else if (line.find("Snake sort iteration") != std::string::npos) {
            snakeCount++;
        }
    }

    std::string bestAlgorithm;
    int minCount = std::min(std::min(insertionCount, countingCount), snakeCount);

    if (minCount == insertionCount) {
        bestAlgorithm = "Insertion sort";
    }
    else if (minCount == countingCount) {
        bestAlgorithm = "Counting sort";
    }
    else {
        bestAlgorithm = "Snake sort";
    }

    std::cout << "Insertion sort iteration = " << insertionCount << std::endl;
    std::cout << "Counting sort iteration = " << countingCount << std::endl;
    std::cout << "Snake sort iteration = " << snakeCount << std::endl;
    std::cout << "The best algorithm: " << bestAlgorithm << std::endl;

    inFile.close();

    return 0;
}
*/













/*
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main()
{
    
    srand(time(NULL));

    
    int lotteryNumbers[7];

    for (int i = 0; i < 7; i++) {
        lotteryNumbers[i] = rand() % 10;
    }

    for (int i = 0; i < 7; i++) {
        cout << lotteryNumbers[i] << " ";
    }

    return 0;
}
*/
























///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//1
/*
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main() {
    string input_str;
    getline(cin, input_str); 

    int min_num = INT_MAX;
    int max_num = INT_MIN;
    int sum = 0;

    stringstream ss(input_str);

    while (!ss.eof()) {
        int num;
        ss >> num;

        if (ss.fail()) {
            ss.clear();
            ss.ignore();
            continue;
        }

        min_num = min(min_num, num);
        max_num = max(max_num, num);
        sum += num;
    }

    cout << "Min: " << min_num << endl;
    cout << "Max: " << max_num << endl;
    cout << "Symma: " << sum << endl;

    return 0;
}*/


//2
/*
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace std;

int main() {
    string input_str;
    getline(cin, input_str); 

    unordered_map<int, int> num_count;

    stringstream ss(input_str);
    int num;
    while (ss >> num) {
        num_count[num]++;
    }

    vector<int> num_list;
    ss.clear();
    ss.str(input_str);
    while (ss >> num) {
        if (num_count[num] > 1) {
            num_list.push_back(num);
        }
    }

    for (int i = 0; i < num_list.size(); i++) {
        cout << num_list[i];
        if (i < num_list.size() - 1) {
            cout << " ";
        }
    }
    cout << endl;

    return 0;
}*/


//3
/*
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

int main() {
    string input_str;
    getline(cin, input_str); 

    vector<int> even_nums; 
    vector<int> odd_nums; 

    stringstream ss(input_str);
    int num;
    char delimiter;
    while (ss >> num) {
        even_nums.push_back(num);
        if (ss >> delimiter) { 
            if (delimiter != ',') {
                cerr << "Ошибка: разделитель не является запятой" << endl;
                return 1;
            }
        }
        else {
            break;
        }
    }

    for (int i = 0; i < even_nums.size(); i++) {
        if (even_nums[i] % 2 == 1) {
            odd_nums.push_back(even_nums[i]);
            even_nums.erase(even_nums.begin() + i); 
            i--;
        }
    }

   
    for (int i = 0; i < odd_nums.size(); i++) {
        cout << odd_nums[i];
        if (i < odd_nums.size() - 1) {
            cout << " ";
        }
    }
    cout << endl;

    for (int i = 0; i < even_nums.size(); i++) {
        cout << even_nums[i];
        if (i < even_nums.size() - 1) {
            cout << " ";
        }
    }
    cout << endl;

    return 0;
}*/


//4
/*
#include <iostream>
#include <string>

using namespace std;

int main() {
    string word;
    cout << "Slovo: ";
    cin >> word;

    char error_letter = '\0';
    int error_count = 0;
    string corrected_word = word;

    for (int i = 1; i < word.length() - 1; i++) {
        if (word[i] == word[i - 1] && word[i] == word[i + 1]) {
            if (error_letter == '\0') {
                error_letter = word[i];
                error_count = 3;
            }
            else if (word[i] == error_letter) {
                error_count += 3;
            }
            else {
                cout << "I cant" << endl;
                return 0;
            }
        }
    }

    if (error_letter != '\0') {
        if (error_count == 3) {
            char correct_letter = (error_letter == 'a') ? 'a' : 'a';
            corrected_word.replace(corrected_word.find(error_letter), 3, 1, correct_letter);
        }
        cout << "Bykva \"" << error_letter << "\" Oshibochno povtor " << error_count / 3 << " raz" << endl;
        cout << "Ispravlennoe: " << corrected_word << endl;
    }
    else {
        cout << "No error" << endl;
    }

    return 0;
}*/


//5
/*
#include <iostream>
#include <string>

int main() {
    std::string input_string;
    std::getline(std::cin, input_string);  
    int index_of_hash = input_string.find('#');  

    if (index_of_hash != std::string::npos) {  
        input_string.erase(index_of_hash, 1); 
    }

    std::cout << "ispravlenno: " << input_string << std::endl;

    std::string longest_word;
    std::string current_word;

    
    for (char c : input_string) {
        if (std::isalpha(c)) {  
            current_word += c;  
        }
        else {
            
            if (current_word.length() > longest_word.length()) {
                longest_word = current_word;  
            }
            current_word = "";  
        }
    }

    
    if (current_word.length() > longest_word.length()) {
        longest_word = current_word;
    }

    std::cout << "dlinnoe: " << longest_word << std::endl;

    return 0;
}
*/


//6
/*
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

int main() {
    std::string input_string;
    std::cout << "Enter a string of words: ";
    std::getline(std::cin, input_string);

    std::vector<std::string> words;
    std::string current_word = "";
    for (char c : input_string) {
        if (isalpha(c)) {
            current_word += c;
        }
        else if (!current_word.empty()) {
            words.push_back(current_word);
            current_word = "";
        }
    }
    if (!current_word.empty()) {
        words.push_back(current_word);
    }

    std::sort(words.begin(), words.end());

    std::cout << "Sorted words:\n";
    for (const std::string& word : words) {
        std::cout << word << '\n';
    }

    return 0;
}
*/

//7
/*
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

bool is_word(const std::string& str) {
    for (char c : str) {
        if (!isalpha(c)) {
            return false;
        }
    }
    return true;
}

bool compare(const std::string& a, const std::string& b) {
    return a < b;
}

int main() {
    std::string input;
    std::getline(std::cin, input);

    std::vector<std::string> words;
    std::string current_word;
    for (char c : input) {
        if (isalpha(c)) {
            current_word += c;
        }
        else {
            if (is_word(current_word)) {
                words.push_back(current_word);
            }
            current_word.clear();
        }
    }
    if (is_word(current_word)) {
        words.push_back(current_word);
    }

    std::vector<std::string> vowels;
    std::vector<std::string> consonants;
    for (const std::string& word : words) {
        if (word.empty()) {
            continue;
        }
        if (word[0] == 'a' || word[0] == 'e' || word[0] == 'i' || word[0] == 'o' || word[0] == 'u' ||
            word[0] == 'A' || word[0] == 'E' || word[0] == 'I' || word[0] == 'O' || word[0] == 'U') {
            vowels.push_back(word);
        }
        else {
            consonants.push_back(word);
        }
    }

    std::sort(vowels.begin(), vowels.end(), compare);
    std::sort(consonants.begin(), consonants.end(), compare);

    for (const std::string& word : vowels) {
        std::cout << word << " ";
    }
    std::cout << std::endl;

    for (const std::string& word : consonants) {
        std::cout << word << " ";
    }
    std::cout << std::endl;

    return 0;
}*/


//8
/*
#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main() {
    vector<string> keys = { "model", "price", "quantity", "size", "color", "discount" };
    vector<string> values = { "XXLDude", "5678.00", "57", "XXL", "black", "12%" };

    for (int i = 0; i < keys.size(); i++) {
        cout << keys[i] << ": " << values[i] << endl;
    }

    return 0;
}
*/
//9
/*
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int main() {
    string input;
    cout << "Enter student name and grades: ";
    getline(cin, input);

    istringstream iss(input);
    string name;
    vector<int> grades;
    iss >> name;
    int grade;
    while (iss >> grade) {
        grades.push_back(grade);
    }

    double sum = 0;
    for (int grade : grades) {
        sum += grade;
    }
    double average = sum / grades.size();

    cout << name << ", average grade: " << average << endl;
    return 0;
}*/


//10
/*
#include <iostream>
#include <string>
#include <vector>

int main() {
    std::string input;
    std::getline(std::cin, input);

    std::vector<std::string> route;
    std::string temp;
    for (int i = 0; i < input.length(); i++) {
        if (input[i] == ' ') {
            route.push_back(temp);
            temp.clear();
        }
        else {
            temp.push_back(input[i]);
        }
    }
    route.push_back(temp);

    int s = 0, n = 0, w = 0, e = 0;
    for (auto i : route) {
        char direction = i[0];
        int distance = std::stoi(i.substr(1));

        if (direction == 's') {
            s += distance;
        }
        else if (direction == 'n') {
            n += distance;
        }
        else if (direction == 'w') {
            w += distance;
        }
        else if (direction == 'e') {
            e += distance;
        }
    }

    std::cout << "x: " << e - w << ", y: " << n - s << std::endl;
    return 0;
}*/



























// 100 б из них счастливый есть есть 50-55 100б рандомно
/*
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

int main()
{
    srand(time(NULL));

    vector<string> tickets;

    for (int i = 0; i < 100; i++) {
        string ticket = "";
        for (int j = 0; j < 7; j++) {
            ticket += to_string(rand() % 10);
        }
        tickets.push_back(ticket);
    }

    ofstream file("C:\\Users\\ASUS\\Desktop\\vvv.txt");
    for (const auto& ticket : tickets) {
        file << ticket << endl;
    }
    file.close();

    vector<string> lucky_tickets;
    for (const auto& ticket : tickets) {
        if (ticket.find("50") != string::npos ||
            ticket.find("51") != string::npos ||
            ticket.find("52") != string::npos ||
            ticket.find("53") != string::npos ||
            ticket.find("54") != string::npos ||
            ticket.find("55") != string::npos) {
            lucky_tickets.push_back(ticket);
        }
    }

    sort(lucky_tickets.rbegin(), lucky_tickets.rend());

    if (!lucky_tickets.empty()) {
        cout << "Lucky ticket: " << lucky_tickets[0] << endl;
    }
    else {
        cout << "No lucky ticket found." << endl;
    }

    return 0;
}*/



//1
/*
4
123456 Ivan
234567 Petro
345678 Anna
456789 Maria
Petro*/
/*
#include <iostream>
#include <string>
#include <map>

using namespace std;

int main() {
    int n;
    cin >> n;

    map<string, string> phonebook;

    for (int i = 0; i < n; i++) {
        string phone, name;
        cin >> phone >> name;
        phonebook[name] = phone;
    }

    string query;
    cin >> query;

    if (phonebook.find(query) != phonebook.end()) {
        cout << phonebook[query] << endl;
    }
    else {
        cout << "No user" << endl;
    }

    return 0;
}*/




//2
/*
5
Ivanov informatics 101
Petrov informatics 102
Sidorov math 103
Kozlov informatics 103
Novikov physics 104
informatics
*/
/*
#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;

int main() {
    int n;
    cin >> n;

    map<string, vector<string>> students;

    for (int i = 0; i < n; i++) {
        string surname, specialty, group;
        cin >> surname >> specialty >> group;
        students[specialty].push_back(surname);
    }

    string query;
    cin >> query;

    if (students.find(query) != students.end()) {
        vector<string> surnames = students[query];
        for (int i = 0; i < surnames.size(); i++) {
            cout << surnames[i];
            if (i < surnames.size() - 1) {
                cout << ", ";
            }
        }
        cout << endl;
    }
    else {
        cout << "No user" << endl;
    }

    return 0;
}*/


//3 
/*
#include <iostream>
#include <string>
#include <map>

using namespace std;

int main() {
    string text;
    getline(cin, text);

    map<char, string> translit = {
        {'А', "A"}, {'Б', "B"}, {'В', "V"}, {'Г', "G"}, {'Д', "D"},
        {'Е', "E"}, {'Ё', "E"}, {'Ж', "ZH"}, {'З', "Z"}, {'И', "I"},
        {'Й', "I"}, {'К', "K"}, {'Л', "L"}, {'М', "M"}, {'Н', "N"},
        {'О', "O"}, {'П', "P"}, {'Р', "R"}, {'С', "S"}, {'Т', "T"},
        {'У', "U"}, {'Ф', "F"}, {'Х', "KH"}, {'Ц', "TC"}, {'Ч', "CH"},
        {'Ш', "SH"}, {'Щ', "SHCH"}, {'Ы', "Y"}, {'Ъ', ""}, {'Ь', ""},
        {'Э', "E"}, {'Ю', "IU"}, {'Я', "IA"},
        {'а', "a"}, {'б', "b"}, {'в', "v"}, {'г', "g"}, {'д', "d"},
        {'е', "e"}, {'ё', "e"}, {'ж', "zh"}, {'з', "z"}, {'и', "i"},
        {'й', "i"}, {'к', "k"}, {'л', "l"}, {'м', "m"}, {'н', "n"},
        {'о', "o"}, {'п', "p"}, {'р', "r"}, {'с', "s"}, {'т', "t"},
        {'у', "u"}, {'ф', "f"}, {'х', "kh"}, {'ц', "tc"}, {'ч', "ch"},
        {'ш', "sh"}, {'щ', "shch"}, {'ы', "y"}, {'ъ', ""}, {'ь', ""},
        {'э', "e"}, {'ю', "iu"}, {'я', "ia"}
    };

    string result;
    for (char c : text) {
        if (translit.find(c) != translit.end()) {
            result += translit[c];
        }
        else {
            result += c;
        }
    }

    cout << result << endl;

    return 0;
}*/






/////////////////////////////////////////////////////////////////////////////////
//Колок 15


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;


struct Respondent {
    int age;
    string gender;
    string education;
    bool answer;
};


void addRespondentToFile(fstream& file, Respondent respondent) {
    file << respondent.age << " " << respondent.gender << " " << respondent.education << " " << respondent.answer << endl;
}


vector<Respondent> readRespondentsFromFile(fstream& file) {
    vector<Respondent> respondents;
    file.clear();
    file.seekg(0, ios::beg);

    int age;
    string gender, education;
    bool answer;
    while (file >> age >> gender >> education >> answer) {
        Respondent respondent = { age, gender, education, answer };
        respondents.push_back(respondent);
    }
    return respondents;
}


void printRespondents(vector<Respondent> respondents) {
    for (auto respondent : respondents) {
        cout << "Age: " << respondent.age << ", Gender: " << respondent.gender << ", Education: " << respondent.education << ", Answer: " << respondent.answer << endl;
    }
}


int countRespondents(const vector<Respondent>& respondents, string gender, int ageThreshold, char ageComparisonOperator, string education, bool answer) {
    int count = 0;
    for (auto respondent : respondents) {
        bool ageCondition = false;
        if (ageComparisonOperator == '>') {
            ageCondition = respondent.age >= ageThreshold;
        }
        else if (ageComparisonOperator == '<') {
            ageCondition = respondent.age <= ageThreshold;
        }
        if (respondent.gender == gender && ageCondition && respondent.education == education && respondent.answer == answer) {
            count++;
        }
    }
    return count;
}



int main() {
   
    fstream file("C:\\Users\\ASUS\\Desktop\\respondents.txt", ios::in | ios::out | ios::app);
    if (!file) {
        cout << "Error opening file" << endl;
        return 1;
    }

    int choice;
    do {
        
        cout << "Menu" << endl;
        cout << "1. Add respondent" << endl;
        cout << "2. Print respondents" << endl;
        cout << "3. Count respondents" << endl;
        cout << "0. Exit" << endl;
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
        case 1: {
           
            Respondent respondent;
            cout << "Enter age: ";
            cin >> respondent.age;
            cout << "Enter gender: ";
            cin >> respondent.gender;
            cout << "Enter education (primary, middle, high): ";
            cin >> respondent.education;
            cout << "Enter answer (0 for NO, 1 for YES): ";
            cin >> respondent.answer;

            
            addRespondentToFile(file, respondent);

            break;
        }
        case 2: {
            
            file.seekg(0, ios::beg);
            vector<Respondent> respondents = readRespondentsFromFile(file);
            printRespondents(respondents);

            break;
        }
        case 3: {
          
            file.seekg(0, ios::beg);
            vector<Respondent> respondents = readRespondentsFromFile(file);

            
            int count;
            cout << "Enter criteria:" << endl;
            cout << "Gender (male/female): ";
            string gender;
            cin >> gender;
            cout << "Age threshold: ";
            int ageThreshold;
            cin >> ageThreshold;
            cout << "Age comparison operator ('>' or '<'): ";
            char ageComparisonOperator;
            cin >> ageComparisonOperator;
            cout << "Education (primary/middle/high): ";
            string education;
            cin >> education;
            cout << "Answer (0 for NO, 1 for YES): ";
            bool answer;
            cin >> answer;

            count = countRespondents(respondents, gender, ageThreshold, ageComparisonOperator, education, answer);

            cout << "Number of respondents: " << count << endl;

            break;
        }

        case 0: {
            
            break;
        }
        default: {
            cout << "Invalid choice" << endl;
        }
        }

    } while (choice != 0);

  
    file.close();

    return 0;
}


