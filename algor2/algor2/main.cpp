//на пр1
//метод деления пополам
/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return x * x * x + x - 1;
}

int main() {
    double a = 0.0; 
    double b = 1.0; 
    double epsilon = 0.01; 
    double x, fx;

    int iteration = 1;

    while ((b - a) > epsilon) {
        x = (a + b) / 2;
        fx = f(x);

        cout << "iter " << iteration << ":" << endl;
        cout << "a = " << a << ", b = " << b << ", x = " << x << ", b-a = " << b - a << endl;
        cout << "F(a) = " << f(a) << ", F(b) = " << f(b) << ", F(x) = " << fx << endl;

        if (fx == 0) {
            break; 
        }
        else if (fx * f(a) < 0) {
            b = x;
        }
        else {
            a = x; 
        }

        iteration++;
    }

    cout << "res koren: " << x << endl;
    return 0;
}*/



// Метод Ньютона (метод касательных)
/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return x * x * x + x - 1;
}

double df(double x) {
    return 3 * x * x + 1;
}

int main() {
    double epsilon = 0.01; 
    double x0 = 1; 

    int iteration = 1;

    while (true) {
        double fx = f(x0);
        double dfx = df(x0);
        double x1 = x0 - fx / dfx;
        double error = abs(x1 - x0);

        cout << "iter " << iteration << ":" << endl;
        cout << "x = " << x0 << ", F(x) = " << fx << ", F'(x) = " << dfx << ", pogreshnost = " << error << endl;

        if (error < epsilon) {
            cout << "res korn: " << x1 << endl;
            break;
        }

        x0 = x1;
        iteration++;
    }

    return 0;
}*/







//ДЗ 1 метод касательных и пополам


/*
#include <iostream>
#include <cmath>

using namespace std;

double equation(double x) {
    return x + log(x);
}

double derivative(double x) {
    return 1.0 + 1.0 / x;
}

// Метод Касательных (метод Ньютона)
double newtonRaphson(double initialGuess, double epsilon, int maxIterations) {
    double x = initialGuess;
    int iterations = 0;

    while (iterations < maxIterations) {
        double f = equation(x);
        double f_prime = derivative(x);

        if (fabs(f_prime) < 1e-10) {
            cout << "BLIZKO K 0. END." << endl;
            break;
        }

        double x_new = x - f / f_prime;
        iterations++;

        double error = fabs(x_new - x); // Погрешность

        cout << "ITER " << iterations << ": x = " << x_new << ", f(x) = " << f << ", POGRESHNOST: " << error << endl;

        if (error < epsilon) {
            cout << "NAIDENO:" << endl;
            cout << "KOREN X: " << x_new << ", ITERS: " << iterations << ", f(x) = " << equation(x_new) << endl;
            return x_new;
        }

        x = x_new;
    }

    cout << "ERROR " << maxIterations << " ." << endl;
    return initialGuess;
}
/*
// Метод половинного деления
double bisection(double a, double b, double epsilon, int maxIterations) {
    int iterations = 0;

    if (equation(a) * equation(b) >= 0) {
        cout << "ERROR." << endl;
        return a;
    }

    while (iterations < maxIterations) {
        double c = (a + b) / 2;
        iterations++;

        double error = fabs(equation(c)); // Погрешность

        cout << "ITERS " << iterations << ": a = " << a << ", b = " << b << ", c = " << c << ", f(c) = " << equation(c) << ", POGRESHNOST: " << error << endl;

        if (error < epsilon) {
            cout << "NAIDENO:" << endl;
            cout << "KOREN X: " << c << ", ITERS: " << iterations << ", f(c) = " << equation(c) << endl;
            return c;
        }

        if (equation(c) * equation(a) < 0) {
            b = c;
        }
        else {
            a = c;
        }
    }

    cout << "ERROR " << maxIterations << " ." << endl;
    return (a + b) / 2;
}

int main() {
    double epsilon = 0.00000001;
    int maxIterations = 100;
    double initialGuess = 0.57;

    cout << "METOD KASATELNIX:" << endl;
    double newtonRoot = newtonRaphson(initialGuess, epsilon, maxIterations);

    cout << endl;

    cout << "METOD DELENIA POPOLAM:" << endl;
    double a = 0.1; // Начальный интервал [a, b]
    double b = 1.0;
    double bisectionRoot = bisection(a, b, epsilon, maxIterations);
    
    return 0;
}
*/



/*
#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return x + log(x);
}

int main() {
    double a = 0.0;
    double b = 1.0;
    double epsilon = 0.01;
    double x, fx;

    int iteration = 1;

    while ((b - a) > epsilon) {
        x = (a + b) / 2;
        fx = f(x);

        cout << "iter " << iteration << ":" << endl;
        cout << "a = " << a << ", b = " << b << ", x = " << x << ", b-a = " << b - a << endl;
        cout << "F(a) = " << f(a) << ", F(b) = " << f(b) << ", F(x) = " << fx << endl;

        if (fx == 0) {
            break;
        }
        else if (fx * f(a) < 0) {
            b = x;
        }
        else {
            a = x;
        }

        iteration++;
    }

    cout << "res koren: " << x << endl;
    return 0;
}*/

/*

#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    return x + log(x);
}

double df(double x) {
    return 1 + 1 / x;
}

int main() {
    double epsilon = 0.01;
    double x0 = 1;

    int iteration = 1;

    while (true) {
        double fx = f(x0);
        double dfx = df(x0);
        double x1 = x0 - fx / dfx;
        double error = abs(x1 - x0);

        cout << "iter " << iteration << ":" << endl;
        cout << "x = " << x0 << ", F(x) = " << fx << ", F'(x) = " << dfx << ", pogreshnost = " << error << endl;

        if (error < epsilon) {
            cout << "res korn: " << x1 << endl;
            break;
        }

        x0 = x1;
        iteration++;
    }

    return 0;
}*/









//
/*
#include <iostream>
#include <cmath>

using namespace std;

const double epsilon = 0.001;

// Функция, заданная уравнением
double f(double x) {
    return 2 * sin(x) - atan(x);
}

int main() {
    double x0 = 2.5; // Начальное приближение на интервале (2.5, 2.6)
    double x1;
    int iteration = 0;

    cout << " x0 = " << x0 << endl;

    do {
        x1 = x0 - f(x0) / (2 * cos(x0) - 1 / (1 + x0 * x0));
        double current_epsilon = fabs(x1 - x0);
        cout << "iter " << iteration << ": x0 = " << x0 << ", x1 = " << x1 << ", epsilon = " << current_epsilon << endl;
        x0 = x1;
        iteration++;
    } while (fabs(x1 - x0) > epsilon);

    cout << "koren: x = " << x1 << endl;

    return 0;
}*/





/////////////////////////////////

//хорд
/*
#include <cmath>
#include <iostream>

double f(double x) {
    return 2 * sin(x) - atan(x);
}

void chordMethod(double a, double b, double eps) {
    double x = a;
    int iteration = 0;
    while (fabs(f(x)) >= eps) {
        std::cout << "iter " << ++iteration << ":\n";
        std::cout << "a = " << a << ", b = " << b << ", x = " << x << ", f(x) = " << f(x) << "\n";
        x = b - (f(b) * (b - a)) / (f(b) - f(a));
        a = b;
        b = x;
    }
    std::cout << "koren x = : " << x << "\n";
}

int main() {
    double a = 2.5;
    double b = 2.6;
    double eps = 0.0001;

    std::cout.precision(5);
    chordMethod(a, b, eps);

    return 0;
}*/



//касательных
/**/
/*
#include <cmath>
#include <iostream>

double f(double x) {
    return 2 * sin(x) - atan(x);
}

double df(double x) {
    return 2 * cos(x) - 1 / (1 + x * x);
}

void tangentMethod(double a, double b, double eps) {
    double x = (a + b) / 2;
    int iteration = 0;
    while (fabs(f(x)) >= eps) {
        std::cout << "iter " << ++iteration << ":\n";
        std::cout << "a = " << a << ", b = " << b << ", x = " << x << ", f(x) = " << f(x) << "\n";
        x = x - f(x) / df(x);
    }
    std::cout << "x= : " << x << "\n";
}

int main() {
    double a = 2.5;
    double b = 2.6;
    double eps = 0.0001;

    std::cout.precision(5);
    tangentMethod(a, b, eps);

    return 0;
}*/




//гауса
/*
#include <iostream>
using namespace std;

int n, i, j, k;
double d, s;

int main()
{
    cout << "Kol-vo yravneniy: " << endl;
    cin >> n;

    double** a = new double* [n];
    for (i = 0; i < n; i++)
        a[i] = new double[n];

    double* b = new double[n];
    double* x = new double[n];

    cout << "Vvod i " << endl;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            cout << "a[" << i + 1 << "," << j + 1 << "]= ";
            cin >> a[i][j];
        }
        cout << "b[" << i + 1 << "]= ";
        cin >> b[i];
    }

    for (k = 0; k < n; k++)
    {
        cout << "Iter " << k + 1 << ":" << endl;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                cout << a[i][j] << "\t";
            }
            cout << "| " << b[i] << endl;
        }
        for (j = k + 1; j < n; j++)
        {
            d = a[j][k] / a[k][k];
            for (i = k; i < n; i++)
            {
                a[j][i] = a[j][i] - d * a[k][i];
            }
            b[j] = b[j] - d * b[k];
        }
    }

    for (k = n - 1; k >= 0; k--)
    {
        x[k] = b[k];
        for (j = k + 1; j < n; j++)
        {
            x[k] -= a[k][j] * x[j];
        }
        x[k] /= a[k][k];
    }

    cout << "Korni sistemy: " << endl;
    for (i = 0; i < n; i++)
        cout << "x[" << i + 1 << "]=" << x[i] << " " << endl;

    for (i = 0; i < n; i++)
    {
        delete[] a[i];
    }
    delete[] a;
    delete[] b;
    delete[] x;

    return 0;
}*/


//крамера
/*
#include <iostream>
using namespace std;

int determinant(int matrix[3][3]);
int determinantX1(int coefMatrix[3][3], int constTermsMatrix[3]);
int determinantX2(int coefMatrix[3][3], int constTermsMatrix[3]);
int determinantX3(int coefMatrix[3][3], int constTermsMatrix[3]);

int main()
{
    int i, j;

    int coefficientsMatrix3x3[3][3];
    int constantTermsMatrix3x1[3];

    cout << "Vvedite koefficienty i sbobodnye chleny " << endl;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            cout << "a[ " << i << "," << j << "]= ";
            cin >> coefficientsMatrix3x3[i][j];
        }
        cout << "b,[ " << i << "]= ";
        cin >> constantTermsMatrix3x1[i];
    }

    int det = determinant(coefficientsMatrix3x3);
    int detX1 = determinantX1(coefficientsMatrix3x3, constantTermsMatrix3x1);
    int detX2 = determinantX2(coefficientsMatrix3x3, constantTermsMatrix3x1);
    int detX3 = determinantX3(coefficientsMatrix3x3, constantTermsMatrix3x1);

    if (det != 0)
    {
        cout << "X1 = " << (float)detX1 / (float)det << endl;
        cout << "X2 = " << (float)detX2 / (float)det << endl;
        cout << "X3 = " << (float)detX3 / (float)det << endl;
    }
    else
        cout << "Sistema ne imejet reshenij " << endl << endl;

    return 0;
}

int determinant(int matrix[3][3])
{
    int a11 = matrix[0][0];
    int a12 = matrix[0][1];
    int a13 = matrix[0][2];
    int a21 = matrix[1][0];
    int a22 = matrix[1][1];
    int a23 = matrix[1][2];
    int a31 = matrix[2][0];
    int a32 = matrix[2][1];
    int a33 = matrix[2][2];

    cout << "Determinant Matrix:" << endl;
    cout << a11 << "\t" << a12 << "\t" << a13 << endl;
    cout << a21 << "\t" << a22 << "\t" << a23 << endl;
    cout << a31 << "\t" << a32 << "\t" << a33 << endl;
    cout << "D= " << (a11 * a22 * a33) + (a12 * a23 * a31) + (a13 * a21 * a32) -
        (a13 * a22 * a31) - (a11 * a23 * a32) - (a12 * a21 * a33) << endl;

    return (a11 * a22 * a33) + (a12 * a23 * a31) + (a13 * a21 * a32) -
        (a13 * a22 * a31) - (a11 * a23 * a32) - (a12 * a21 * a33);
}

int determinantX1(int coefMatrix[3][3], int constTermsMatrix[3])
{
    int a12 = coefMatrix[0][1];
    int a13 = coefMatrix[0][2];
    int a22 = coefMatrix[1][1];
    int a23 = coefMatrix[1][2];
    int a32 = coefMatrix[2][1];
    int a33 = coefMatrix[2][2];
    int c1 = constTermsMatrix[0];
    int c2 = constTermsMatrix[1];
    int c3 = constTermsMatrix[2];

    cout << "DeterminantX1 Matrix:" << endl;
    cout << a12 << "\t" << a13 << "\t" << c1 << endl;
    cout << a22 << "\t" << a23 << "\t" << c2 << endl;
    cout << a32 << "\t" << a33 << "\t" << c3 << endl;
    cout << "D= " << (c1 * a22 * a33) + (a12 * a23 * c3) + (a13 * c2 * a32) -
        (a13 * a22 * c3) - (c1 * a23 * a32) - (a12 * c2 * a33) << endl;

    return (c1 * a22 * a33) + (a12 * a23 * c3) + (a13 * c2 * a32) -
        (a13 * a22 * c3) - (c1 * a23 * a32) - (a12 * c2 * a33);
}

int determinantX2(int coefMatrix[3][3], int constTermsMatrix[3])
{
    int a11 = coefMatrix[0][0];
    int a13 = coefMatrix[0][2];
    int a21 = coefMatrix[1][0];
    int a23 = coefMatrix[1][2];
    int a31 = coefMatrix[2][0];
    int a33 = coefMatrix[2][2];
    int c1 = constTermsMatrix[0];
    int c2 = constTermsMatrix[1];
    int c3 = constTermsMatrix[2];

    cout << "DeterminantX2 Matrix:" << endl;
    cout << a11 << "\t" << a13 << "\t" << c1 << endl;
    cout << a21 << "\t" << a23 << "\t" << c2 << endl;
    cout << a31 << "\t" << a33 << "\t" << c3 << endl;
    cout << "D= " << (a11 * c2 * a33) + (c1 * a23 * a31) + (a13 * a21 * c3) -
        (a13 * c2 * a31) - (a11 * a23 * c3) - (c1 * a21 * a33) << endl;

    return (a11 * c2 * a33) + (c1 * a23 * a31) + (a13 * a21 * c3) -
        (a13 * c2 * a31) - (a11 * a23 * c3) - (c1 * a21 * a33);
}

int determinantX3(int coefMatrix[3][3], int constTermsMatrix[3])
{
    int a11 = coefMatrix[0][0];
    int a12 = coefMatrix[0][1];
    int a21 = coefMatrix[1][0];
    int a22 = coefMatrix[1][1];
    int a31 = coefMatrix[2][0];
    int a32 = coefMatrix[2][1];
    int c1 = constTermsMatrix[0];
    int c2 = constTermsMatrix[1];
    int c3 = constTermsMatrix[2];

    cout << "DeterminantX3 Matrix:" << endl;
    cout << a11 << "\t" << a12 << "\t" << c1 << endl;
    cout << a21 << "\t" << a22 << "\t" << c2 << endl;
    cout << a31 << "\t" << a32 << "\t" << c3 << endl;
    cout << "D= " << (a11 * a22 * c3) + (a12 * c2 * a31) + (c1 * a21 * a32) -
        (c1 * a22 * a31) - (a11 * c2 * a32) - (a12 * a21 * c3) << endl;

    return (a11 * a22 * c3) + (a12 * c2 * a31) + (c1 * a21 * a32) -
        (c1 * a22 * a31) - (a11 * c2 * a32) - (a12 * a21 * c3);
}
*/





//////////////////////////////////////////
//////////////////////////////////////////
/////////////////////////////////////////










/*
15

x1+ 0,5x2 =3
2x1-5x2+ x3=1
    x2+ 8x3- 2x4=5
        1.5x3 -6x4=4
*/

//3.3. Метод прогонки

/*
#include <iostream>
#include <vector>
using namespace std;

// Функция для решения системы вида Ax = d, где A - трехдиагональная матрица
vector<double> solve_tridiagonal(vector<double> a, vector<double> b, vector<double> c, vector<double> d) {
    int n = a.size();
    vector<double> x(n);
    vector<double> u(n); 
    vector<double> v(n);

    // Прямой ход прогонки
    u[0] = -c[0] / b[0];
    v[0] = d[0] / b[0];
    cout << "U[" << 1 << "] = " << u[0] << ", V[" << 1 << "] = " << v[0] << endl;
    for (int i = 1; i < n; i++) {
        double denom = a[i] * u[i - 1] + b[i];
        u[i] = -c[i] / denom;
        v[i] = (d[i] - a[i] * v[i - 1]) / denom;
        cout << "U[" << i + 1 << "] = " << u[i] << ", V[" << i + 1 << "] = " << v[i] << endl;
    }

    // Обратный ход прогонки
    x[n - 1] = v[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = u[i] * x[i + 1] + v[i];
    }

    return x;
}


int main() {
    
    vector<double> a = { 0, 2, 0, 1.5 }; 
    vector<double> b = { 1, -5, 8, -6 };   
    vector<double> c = { 0.5, 1, -2, 0 };  
    vector<double> d = { 3, 1, 5, 4 }; 
    cout << "a: ";
    for (double ai : a) {
        cout << ai << " ";
    }
    cout << endl;

    cout << "b: ";
    for (double bi : b) {
        cout << bi << " ";
    }
    cout << endl;

    cout << "c: ";
    for (double ci : c) {
        cout << ci << " ";
    }
    cout << endl;

    cout << "d: ";
    for (double di : d) {
        cout << di << " ";
    }
    cout << endl;

    
    vector<double> x = solve_tridiagonal(a, b, c, d);

  
    cout << "Res:" << endl;
    for (int i = 0; i < x.size(); i++) {
        cout << "x[" << i + 1 << "] = " << x[i] << endl;
    }

   
    vector<double> true_solution = { 2.552536231884058, 0.8949275362318841, 0.3695652173913043, -0.5742753623188406 };

    
    cout << "Pogreshnost:" << endl;
    for (int i = 0; i < x.size(); i++) {
        double error = abs(x[i] - true_solution[i]);
        cout << "x[" << i + 1 << "] = " << x[i] << ", Pogreshnost = " << error << endl;
    }


    return 0;
}*/




//метод Якоби
/*
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Функция для вычисления новых значений переменных на следующей итерации
std::vector<double> calculateNextIteration(const std::vector<double>& x) {
    std::vector<double> newX(x.size());
    newX[0] = (3 - 0.5 * x[1]) / 1.0;
    newX[1] = (1 - 2 * x[0] + x[2]) / -5.0;
    newX[2] = (5 - x[1] + 2 * x[3]) / 8.0;
    newX[3] = (4 - 1.5 * x[2]) / -6.0;
    return newX;
}

// Функция для проверки сходимости метода
bool isConverged(const std::vector<double>& x, const std::vector<double>& newX, double tolerance) {
    for (size_t i = 0; i < x.size(); ++i) {
        if (std::abs(x[i] - newX[i]) > tolerance) {
            return false;
        }
    }
    return true;
}

// Функция для вывода вектора на экран
void print_vector(const vector<double>& x, const string& name) {
    int n = x.size(); 
    cout << name << " = (";
    for (int i = 0; i < n; i++) {
        cout << x[i];
        if (i < n - 1) cout << ", ";
    }
    cout << ")" << endl;
}

int main() {
    const double tolerance = 0.0000001;
    std::vector<double> x = { 0.0, 0.0, 0.0, 0.0 };
    std::vector<double> newX;

    int maxIterations = 1000;

    for (int iteration = 1; iteration <= maxIterations; ++iteration) {
        newX = calculateNextIteration(x);

        
        cout << "iter " << iteration << ":" << endl;
        print_vector(newX, "x");

        if (isConverged(x, newX, tolerance)) {
            std::cout << "kol-vo " << iteration << std::endl;
            break;
        }

        x = newX;
    }

    std::vector<double> true_solution = { 2.552536231884058, 0.8949275362318841, 0.3695652173913043, -0.5742753623188406 };

    cout << "Pogreshnost:" << endl;
    for (size_t i = 0; i < x.size(); ++i) {
        double relative_error = std::abs(x[i] - true_solution[i]);
        cout << "x[" << i + 1 << "] = " << x[i] << ", Pogreshnost = " << relative_error << endl;
    }

    return 0;
}*/


//Зейделя

/*
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

const double EPSILON = 0.0000001;
const int MAX_ITERATIONS = 1000;

// Функция для вычисления новых значений переменных методом Зейделя
void gaussSeidel(vector<vector<double>>& A, vector<double>& b, vector<double>& x) {
    int n = x.size();
    vector<double> newX(n);

    for (int iteration = 0; iteration < MAX_ITERATIONS; ++iteration) {
        for (int i = 0; i < n; ++i) {
            double sum1 = 0.0;
            double sum2 = 0.0;

            for (int j = 0; j < n; ++j) {
                if (j < i) {
                    sum1 += A[i][j] * newX[j];
                }
                else if (j > i) {
                    sum2 += A[i][j] * x[j];
                }
            }

            newX[i] = (b[i] - sum1 - sum2) / A[i][i];
        }

        // Проверяем критерий останова
        double maxDiff = 0.0;
        for (int i = 0; i < n; ++i) {
            maxDiff = max(maxDiff, abs(newX[i] - x[i]));
            x[i] = newX[i];
        }

        cout << "iters " << iteration + 1 << ":" << endl;
        for (int i = 0; i < x.size(); ++i) {
            cout << "x" << i + 1 << " = " << x[i] << endl;
        }

        if (maxDiff < EPSILON) {
            break;
        }
    }
}

int main() {
    vector<vector<double>> A = { {1, 0.5, 0, 0},
                                {2, -5, 1, 0},
                                {0, 1, 8, -2},
                                {0, 0, 1.5, -6} };

    vector<double> b = { 3, 1, 5, 4 };
    vector<double> x = { 0, 0, 0, 0 }; 

    gaussSeidel(A, b, x);

    cout << "Otvet:" << endl;
    for (int i = 0; i < x.size(); ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }

    vector<double> true_solution = { 2.552536231884058, 0.8949275362318841, 0.3695652173913043, -0.5742753623188406 };
    cout << "Pogreshnost:" << endl;
    for (int i = 0; i < x.size(); ++i) {
        double error = abs(x[i] - true_solution[i]);
        cout << "x" << i + 1 << " = " << x[i] << ", Pogreshnost = " << error << endl;
    }

    return 0;
}
*/








//Задача мин стоим 2 задача минимальная сумм





/*
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>

using namespace std;

const int numRows = 3;
const int numCols = 4;

// Заполняет матрицу перевозок на основе исходных данных
void fillTransportMatrix(vector<vector<int>>& transportMatrix) {
    transportMatrix = {
        {7, 8, 1, 2},
        {3, 4, 5, 7},
        {9, 9, 8, 6}
    };
}

// Выводит матрицу или таблицу в консоль
void printMatrix(const vector<vector<int>>& matrix) {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            cout << setw(4) << matrix[i][j] << '\t';
        }
        cout << endl;
    }
    cout << endl;
}

// Находит оптимальное решение задачи
int solveTransportProblem(vector<vector<int>>& transportMatrix) {
    vector<int> supply = { 240, 330, 150 };
    vector<int> demand = { 230, 190, 240, 60 };

    int totalCost = 0;

    while (true) {
        cout << "Iter:" << endl;
        printMatrix(transportMatrix);

        int minCost = INT_MAX;
        int minI = -1, minJ = -1;

        // Находим ячейку с минимальной стоимостью перевозки
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                if (transportMatrix[i][j] != 0) {
                    if (transportMatrix[i][j] < minCost) {
                        minCost = transportMatrix[i][j];
                        minI = i;
                        minJ = j;
                    }
                }
            }
        }

        // Если не нашли, то оптимальное решение найдено
        if (minI == -1) {
            break;
        }

        // Находим количество товаров, которое можно перевезти в этой ячейке
        int minAmount = min(supply[minI], demand[minJ]);

        // Обновляем остатки предложения и спроса
        supply[minI] -= minAmount;
        demand[minJ] -= minAmount;

        // Заполняем ячейку в решении и увеличиваем общую стоимость
        transportMatrix[minI][minJ] = minAmount;
        totalCost += minCost * minAmount;
    }


    return totalCost;
}

int main() {
    vector<vector<int>> transportMatrix;
    fillTransportMatrix(transportMatrix);

    cout << "Do resh:" << endl;
    printMatrix(transportMatrix);

    int minCost = solveTransportProblem(transportMatrix);

    cout << "Min zatrat: " << minCost << endl;

    return 0;
}*/





//задача 1
/*
#include <iostream>
#include <vector>

using namespace std;

const int n = 4;  // число неизвестных
const int m = 3;  // число ограничений
const int m1 = 0; // последняя строка равенств
const int m2 = 1; // последняя строка неравенств вида >=

int Bi[50]; // искусственный базис

vector<vector<double>> a(m, vector<double>(50));
vector<double> b(m);
vector<double> c(50);
vector<double> x(50);
vector<double> e(50);
vector<double> cb(50);

void SimplexMethod() {
    int m21 = m2 - m1 + n;
    int nm1 = n + m - m1;
    int n1 = n + m - m1 + m2;
    double mb = 12345.0;

    // Создание искусственного базиса
    for (int i = 0; i < m2; ++i) {
        cb[i] = mb;
        Bi[i] = nm1 + i;
    }

    for (int i = m2; i < m; ++i) {
        Bi[i] = m21 + i - m2;
        cb[i] = 0.0;
    }

    for (int i = 0; i < n1; ++i) {
        x[i] = 0.0;
    }

    cout << "Rechenie zadachi:" << endl;

    // Применяем симплексный метод, вычисляем оценки
    while (true) {
        for (int j = 0; j < n1; ++j) {
            double s0 = 0.0;
            for (int i = 0; i < m; ++i) {
                s0 += cb[i] * a[i][j];
            }
            e[j] = s0 - c[j];
        }

        double max = e[0];
        int j0 = 0;
        for (int i = 1; i < n1; ++i) {
            if (e[i] > max) {
                max = e[i];
                j0 = i;
            }
        }

        // Получили столбец с максимальной оценкой
        if (max <= 0) {
            for (int i = 0; i < m; ++i) {
                x[Bi[i]] = b[i];
            }
            break;
        }

        double s1 = a[0][j0];
        for (int i = 1; i < m; ++i) {
            if (s1 < a[i][j0]) {
                s1 = a[i][j0];
            }
        }

        if (s1 <= 0) {
            cout << "Net optimalnogo plana, f neogranichena" << endl;
            break;
        }

        s1 = mb;
        int i0 = 0;
        for (int i = 0; i < m; ++i) {
            if (a[i][j0] > 0) {
                double s = b[i] / a[i][j0];
                if (s < s1) {
                    s1 = s;
                    i0 = i;
                }
            }
        }

        // Главный элемент a[i0][j0]
        double s0 = a[i0][j0];
        Bi[i0] = j0;
        for (int j = 0; j < n1; ++j) {
            a[i0][j] = a[i0][j] / s0;
        }
        b[i0] = b[i0] / s0;

        for (int i = 0; i < m; ++i) {
            if (i != i0) {
                s1 = -a[i][j0];
                b[i] = b[i] + b[i0] * s1;
                for (int j = 0; j < n1; ++j) {
                    a[i][j] = a[i][j] + a[i0][j] * s1;
                }
            }
        }

        cb[i0] = c[j0];
    }
}

int main() {
    

    cout << "Vvedite matrisy A " << m << "x" << n << " postrochno:" << endl;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> a[i][j];
        }
    }

    cout << "Vvedite v vide stroki vector b, sost iz " << m << " komponentov:" << endl;
    for (int i = 0; i < m; ++i) {
        cin >> b[i];
    }

    cout << "Vvedite vector c, sost iz " << n << " komponentov:" << endl;
    for (int i = 0; i < n; ++i) {
        cin >> c[i];
    }

    SimplexMethod();

    cout << "Res dla vypyska po resyrs:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x[" << i + 1 << "] = " << x[i] << endl;
    }
    cout << "f = " << 65 * x[0] + 70 * x[1] + 60 * x[2] + 120 * x[3];
    
    return 0;
}*/

















/////////////////////////////
//Метод прямоугольников

/*
#include <iostream>
#include <cmath>

using namespace std;

// Функция, которую нужно проинтегрировать
double f(double x) {
    return 3 * x * x - x - 1;
}

// Функция для вычисления интеграла методом прямоугольников
double rectangularIntegration(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        double x_i = a + i * h + h / 2; // Средняя точка на подинтервале
        double f_x_i = f(x_i); // Значение функции в средней точке
        sum += f_x_i;

        cout << "Iteration " << i + 1 << ": x = " << x_i << ", f(x) = " << f_x_i << ", c = " << sum << ", h = " << h << ", I = " << h * sum << endl;
    }

    return h * sum;
}

int main() {
    double a = 1.0;  // Нижний предел интегрирования
    double b = 3.0;  // Верхний предел интегрирования
    int n = 8;       // Количество подинтервалов

    double result = rectangularIntegration(a, b, n);
    cout << "Res: " << result << endl;

    double trueValue = 20.0; // Истинное значение интеграла

    double error = abs(result - trueValue); // Погрешность

    cout << "True Value: " << trueValue << endl;
    cout << "Pogreshnost: " << error << endl;

    return 0;
}*/







//Метод трапеций
/*
#include <iostream>
#include <cmath>

using namespace std;

// Функция, которую нужно проинтегрировать
double f(double x) {
    return 3 * x * x - x - 1;
}

// Функция для вычисления интеграла методом трапеций
double trapezoidalIntegration(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;

    sum += f(a) / 2.0; // Половина высоты первой трапеции

    for (int i = 1; i < n; i++) {
        double x_i = a + i * h;
        double f_x_i = f(x_i); // Значение функции в точке x_i
        sum += f_x_i;

        cout << "Iter " << i << ": x = " << x_i << ", f(x) = " << f_x_i << ", cf = " << sum << ", h = " << h << ", I = " << h * sum << endl;
    }

    sum += f(b) / 2.0; // Половина высоты последней трапеции

    return h * sum;
}

int main() {
    double a = 1.0;  // Нижний предел интегрирования
    double b = 3.0;  // Верхний предел интегрирования
    int n = 8;     // Количество трапеций

    double result = trapezoidalIntegration(a, b, n);
    cout << "Res: " << result << endl;

    double trueValue = 20.0; // Истинное значение интеграла

    double error = abs(result - trueValue); // Погрешность

    cout << "True Value: " << trueValue << endl;
    cout << "Pogreshnost: " << error << endl;

    return 0;
}*/





//Метод Симпсона
/*

#include <iostream>
#include <cmath>

using namespace std;

// Функция, которую нужно проинтегрировать
double f(double x) {
    return 3 * x * x - x - 1;
}

// Функция для вычисления интеграла методом Симпсона
double simpsonIntegration(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    double intermediateResult = sum;

    for (int i = 1; i < n; i += 2) {
        double x_i = a + i * h;
        intermediateResult += 4 * f(x_i);
    }

    for (int i = 2; i < n - 1; i += 2) {
        double x_i = a + i * h;
        intermediateResult += 2 * f(x_i);
    }

    intermediateResult *= (h / 3);
    sum = intermediateResult;

    // Вывод промежуточных результатов
    for (int i = 0; i < n; i++) {
        double x_i = a + i * h;
        double f_x_i = f(x_i);
        cout << "Iter " << i << ": x = " << x_i << ", f(x) = " << f_x_i
            << ", cf = " << intermediateResult
            << ", h = " << h
            << ", I = " << sum << endl;
    }

    return sum;
}

int main() {
    double a = 1.0;  // Нижний предел интегрирования
    double b = 3.0;  // Верхний предел интегрирования
    int n = 40;       // Количество подинтервалов (чётное число)

    double result = simpsonIntegration(a, b, n);
    cout << "Res: " << result << endl;

    double trueValue = 20.0; // Истинное значение интеграла

    double error = abs(result - trueValue); // Погрешность

    cout << "True Value: " << trueValue << endl;
    cout << "Pogreshnost: " << error << endl;

    return 0;
}*/









//Лагранжа
/*

#include <iostream>
#include <vector>

using namespace std;

class LagrangePolynomial {
public:
    double calculate(double x) const {
        double result = 0.0;
        cout << "L(" << x << ") = ";
        for (size_t i = 0; i < x_values.size(); ++i) {
            double term = y_values[i];
            for (size_t j = 0; j < x_values.size(); ++j) {
                if (j != i) {
                    term *= (x - x_values[j]) / (x_values[i] - x_values[j]);
                }
            }
            result += term;

           
            cout << term;
            if (i < x_values.size() - 1) {
                cout << " + ";
            }
        }
        cout << endl;
        return result;
    }

    void buildEquation() const {
        cout << "L(x) = ";
        for (size_t i = 0; i < x_values.size(); ++i) {
            cout << y_values[i];
            for (size_t j = 0; j < x_values.size(); ++j) {
                if (j != i) {
                    cout << " * (x - " << x_values[j] << ") / (" << x_values[i] - x_values[j] << ")";
                }
            }
            if (i < x_values.size() - 1) {
                cout << " + ";
            }
        }
        cout << endl;
    }

    void addPoint(double x, double y) {
        x_values.push_back(x);
        y_values.push_back(y);
    }

private:
    vector<double> x_values;
    vector<double> y_values;
};

int main() {
    LagrangePolynomial polynomial;

    polynomial.addPoint(1, 1);
    polynomial.addPoint(2, 0);
    polynomial.addPoint(3, 1);

    
    polynomial.calculate(0);


    polynomial.buildEquation();

    return 0;
}

*/




//Ньютона
/*
#include <iostream>
#include <vector>

using namespace std;

class NewtonPolynomial {
public:
    
    double dividedDifference(size_t i, size_t j) const {
        if (i == j) {
            return y_values[i];
        }
        else {
            return (dividedDifference(i + 1, j) - dividedDifference(i, j - 1)) / (x_values[j] - x_values[i]);
        }
    }

    
    double calculate(double x) const {
        double result = 0.0;
        cout << "P(" << x << ") = ";
        for (size_t i = 0; i < x_values.size(); ++i) {
            double term = dividedDifference(0, i);
            for (size_t j = 0; j < i; ++j) {
                term *= (x - x_values[j]);
            }
            result += term;

           
            cout << term;
            if (i < x_values.size() - 1) {
                cout << " + ";
            }
        }
        cout << endl;
        return result;
    }

    
    void buildEquation() const {
        cout << "P(x) = ";
        for (size_t i = 0; i < x_values.size(); ++i) {
            cout << dividedDifference(0, i);
            for (size_t j = 0; j < i; ++j) {
                cout << " * (x - " << x_values[j] << ")";
            }
            if (i < x_values.size() - 1) {
                cout << " + ";
            }
        }
        cout << endl;
    }

   
    void addPoint(double x, double y) {
        x_values.push_back(x);
        y_values.push_back(y);
    }

private:
    vector<double> x_values;
    vector<double> y_values;
};

int main() {
    NewtonPolynomial polynomial;

  
    polynomial.addPoint(1, 1);
    polynomial.addPoint(2, 0);
    polynomial.addPoint(3, 1);

 
    polynomial.calculate(0);

 
    polynomial.buildEquation();

    return 0;
}
*/







//Эйлера


#include <iostream>
#include <cmath>

double dydx(double x, double y) {
    return 5 + x - y;
}

void eulerMethod(double x0, double y0, double x_end, double h) {
    double x = x0;
    double y = y0;

    std::cout << "x\ty" << std::endl;
    std::cout << x << "\t" << y << std::endl;

    while (x < x_end) {
        y += h * dydx(x, y);
        x += h;
        std::cout << x << "\t" << y << std::endl;
    }
}

int main() {
    double x0 = 2.0;  
    double y0 = 1.0;  
    double x_end = 4.0;  
    double h = 0.5;  

    eulerMethod(x0, y0, x_end, h);

    return 0;
}



//Эйлера mod
/*
#include <iostream>
#include <cmath>

double dydx(double x, double y) {
    return 5 + x - y;
}

void modifiedEulerMethod(double x0, double y0, double x_end, double h) {
    double x = x0;
    double y = y0;

    std::cout << "x\ty" << std::endl;
    std::cout << x << "\t" << y << std::endl;

    while (x < x_end) {
        double k1 = h * dydx(x, y);
        double k2 = h * dydx(x + h, y + k1);

        y += 0.5 * (k1 + k2);
        x += h;

        std::cout << x << "\t" << y << std::endl;
    }
}

int main() {
    double x0 = 2.0;  
    double y0 = 1.0;  
    double x_end = 4.0;  
    double h = 0.5;  

    modifiedEulerMethod(x0, y0, x_end, h);

    return 0;
}
*/





//м Рунге-Кутта
/*
#include <iostream>
#include <cmath>

double dydx(double x, double y) {
    return 5 + x - y;
}

void rungeKuttaMethod(double x0, double y0, double x_end, double h) {
    double x = x0;
    double y = y0;

    std::cout << "x\ty" << std::endl;
    std::cout << x << "\t" << y << std::endl;

    while (x < x_end) {
        double k1 = h * dydx(x, y);
        double k2 = h * dydx(x + 0.5 * h, y + 0.5 * k1);
        double k3 = h * dydx(x + 0.5 * h, y + 0.5 * k2);
        double k4 = h * dydx(x + h, y + k3);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        x += h;

        std::cout << x << "\t" << y << std::endl;
    }
}

int main() {
    double x0 = 2.0;  // начальное значение x
    double y0 = 1.0;  // начальное значение y при x = 2
    double x_end = 4.0;  // конечное значение x
    double h = 0.5;  // шаг

    rungeKuttaMethod(x0, y0, x_end, h);

    return 0;
}
*/










//stat   на пайтоне



import numpy as np
from scipy.stats import norm
from scipy.stats import chi2

# Ваши данные
data = np.array([2.316, 2.249, 1.829, 1.547, 2.294, 2.151, 1.066, 2.793, 2.716, 3.044, 2.811, 2.617, 1.755, 1.914, 2.685, 2.791, 2.358, 1.789, 1.844, 2.136, 2.864, 3.958, 2.269, 2.565, 2.557, 3.891, 3.479, -0.294, 2.662, 1.559, 1.526, 2.089, 1.301, 2.347, 1.125, 3.9, -0.726, 2.87, 4.179, 3.085, 3.176, 1.029, 0.424, 2.464, 3.474, 3.165, 1.237, 2.509, 2.496, 1.198])

# 1) Объем выборки
sample_size = len(data)
print(f"1) Объем выборки: {sample_size}")

# 2) Максимальное и минимальное значение
max_value = np.max(data)
min_value = np.min(data)
print(f"2) Максимальное значение: {max_value}, Минимальное значение: {min_value}")

# 3) Размах
range_value = max_value - min_value
print(f"3) Размах: {range_value}")

# 4) Количество интервалов в гистограмме по правилу Стурсерса
num_intervals = int(np.ceil(1 + 3.322 * np.log10(sample_size)))
print(f"4) Количество интервалов в гистограмме по правилу Стурсерса: {num_intervals}")

# 5) Определение длины интервала
interval_length = range_value / num_intervals
print(f"5) Длина интервала: {interval_length}")

# 6) Интервалы, частоты, относительные частоты, высоты столбцов, средние каждого интервала
hist, bin_edges = np.histogram(data, bins = num_intervals)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
relative_freq = hist / sample_size
height = relative_freq / interval_length

for i in range(num_intervals) :
    print(f"6) Интервал {i+1}: {bin_edges[i]:.3f} - {bin_edges[i+1]:.3f}, Частота: {hist[i]}, Относительная частота: {relative_freq[i]:.3f}, Высота столбца: {height[i]:.3f}, Среднее: {bin_centers[i]:.3f}")

    # 7) Вычисление выборочного среднего
    sample_mean = np.mean(data)
    print(f"7) Выборочное среднее: {sample_mean}")

    # 8) Вычисление выборочной дисперсии
    sample_variance = np.var(data, ddof = 0)
    print(f"8) Выборочная дисперсия: {sample_variance}")

    # 9) Вычисление исправленной дисперсии
    corrected_variance = np.var(data, ddof = 1)
    print(f"9) Исправленная дисперсия: {corrected_variance}")

    # 10) Вычисление исправленного среднего квадратического отклонения
    std_deviation = np.std(data, ddof = 1)
    print(f"10) Исправленное среднее квадратическое отклонение: {std_deviation}")

    # 11) Функция плотности нормального распределения
    normal_pdf = norm.pdf(bin_centers, loc = sample_mean, scale = std_deviation)
    print(f"11) Функция плотности нормального распределения: {normal_pdf}")

    # 12) Уровень значимости
    alpha = 0.05
    print(f"12) Уровень значимости: {alpha}")

    # 13) Статистика хи - квадрат Пирсона
    expected_freq = sample_size * normal_pdf * interval_length
    chi_square_stat = np.sum((hist - expected_freq) * *2 / expected_freq)
    chi_square_critical = chi2.ppf(1 - alpha, df = num_intervals - 1)
    print(f"13) Статистика хи-квадрат Пирсона: {chi_square_stat}")
    print(f"   Критическое значение: {chi_square_critical}")

    # 14) Заданная надежность
    confidence_level = 1 - alpha
    print(f"14) Заданная надежность: {confidence_level}")

    # 15) Доверительные границы для математического ожидания
    conf_interval_mean = norm.interval(confidence_level, loc = sample_mean, scale = std_deviation / np.sqrt(sample_size))
    print(f"15) Доверительные границы математического ожидания: {conf_interval_mean}")

    # 16) Доверительные границы для оценки дисперсии
    conf_interval_var = ((sample_size - 1) * corrected_variance / chi2.ppf(1 - alpha / 2, df = sample_size - 1),
        (sample_size - 1) * corrected_variance / chi2.ppf(alpha / 2, df = sample_size - 1))
    print(f"16) Доверительные границы оценки дисперсии: {conf_interval_var}")













// квадрат1
/*
#include <iostream>
#include <vector>
#include <cmath>

struct DataPoint {
    double ti, yi;
};

int main() {
    // Заданные значения коэффициентов
    double k = -9.86;
    double b = 167.85;

    // Заданные значения ti и yi
    std::vector<DataPoint> dataPoints = {
        {10, 68},
        {10.6, 64},
        {11, 59},
        {12, 52},
        {12.5, 45},
        {12.8, 42},
        {13, 38},
        {13.2, 37},
        {13.3, 35},
        {13.7, 34}
    };

    // Инициализация переменных для суммирования
    double sumTi = 0, sumYi = 0, sumTiSquared = 0, sumTiYi = 0, sumDiSquared = 0;

    // Вывод заголовка
    std::cout << "ti\t\tyi\t\tti^2\t\tti*yi\t\ty,ras\t\tdi\t\t\di^2\n";

    // Вычисление и вывод значений для каждой точки
    for (const auto& point : dataPoints) {
        double ti = point.ti;
        double yi = point.yi;

        double predictedY = k * ti + b;
        double di = yi - predictedY;
        double diSquared = di * di;

        // Вывод значений для текущей точки
        std::cout << ti << "\t\t" << yi << "\t\t" << ti * ti << "\t\t" << ti * yi << "\t\t" << predictedY << "\t\t" << di << "\t\t" << diSquared << "\n";

        // Суммирование значений
        sumTi += ti;
        sumYi += yi;
        sumTiSquared += ti * ti;
        sumTiYi += ti * yi;
        sumDiSquared += diSquared;
    }

    // Вывод суммированных значений
    std::cout << "Sum of ti: " << sumTi << "\n";
    std::cout << "Sum of yi: " << sumYi << "\n";
    std::cout << "Sum of ti^2: " << sumTiSquared << "\n";
    std::cout << "Sum of ti*yi: " << sumTiYi << "\n";
    std::cout << "Sum of di^2: " << sumDiSquared << "\n";
    std::cout << "k: " << k << "\n";
    std::cout << "b: " << b << "\n";

    return 0;
}
*/





///квадрат2
/*
#include <iostream>
#include <vector>
#include <iomanip> // Для использования std::fixed

struct DataPoint {
    double xi, yi;
};

int main() {
    // Заданные значения коэффициентов
    double b = 83.73561744;
    double k = -0.902161964;
    double a0 = 72.85684477;
    double a1 = 0.116231688;
    double a2 = -0.022453616;

    // Заданные значения xi и yi
    std::vector<DataPoint> dataPoints = {
        {14, 70},
        {15, 69.5},
        {17, 68.5},
        {19, 67.5},
        {20, 66},
        {21, 65},
        {22, 64.5},
        {24, 62.5},
        {27, 60},
        {32, 53.5}
    };

    // Инициализация переменных для суммирования
    double sumXi = 0, sumYi = 0, sumXiSquared = 0, sumXiCubed = 0, sumXiYi = 0, sumXiSquaredYi = 0, sumYRas = 0;
    double sumXiToFourth = 0;

    // Вывод заголовка
    std::cout << "xi\t\tyi\t\txi^2\t\t\txi^3\t\t\txi^4\t\t\txi*yi\t\t\txi^2*yi\t\t\ty,ras\n";

    // Вычисление и вывод значений для каждой точки
    for (const auto& point : dataPoints) {
        double xi = point.xi;
        double yi = point.yi;

        double xiSquared = xi * xi;
        double xiCubed = xi * xi * xi;
        double xiToFourth = xi * xi * xi * xi;
        double xiYi = xi * yi;
        double xiSquaredYi = xiSquared * yi;
        double yRas = a0 + a1 * xi + a2 * xiSquared;

        // Вывод значений для текущей точки
        std::cout << std::fixed << std::setprecision(4) << xi << "\t\t" << yi << "\t\t" << xiSquared << "\t\t" << xiCubed << "\t\t" << xiToFourth << "\t\t" << xiYi << "\t\t" << xiSquaredYi << "\t\t" << yRas << "\n";

        // Суммирование значений
        sumXi += xi;
        sumYi += yi;
        sumXiSquared += xiSquared;
        sumXiCubed += xiCubed;
        sumXiToFourth += xiToFourth;
        sumXiYi += xiYi;
        sumXiSquaredYi += xiSquaredYi;
        sumYRas += yRas;
    }

    // Вывод суммированных значений
    std::cout << "\nSum of xi: " << sumXi << "\n";
    std::cout << "Sum of yi: " << sumYi << "\n";
    std::cout << "Sum of xi^2: " << sumXiSquared << "\n";
    std::cout << "Sum of xi^3: " << sumXiCubed << "\n";
    std::cout << "Sum of xi^4: " << sumXiToFourth << "\n";
    std::cout << "Sum of xi*yi: " << sumXiYi << "\n";
    std::cout << "Sum of xi^2*yi: " << sumXiSquaredYi << "\n";
    std::cout << "Sum of y,ras: " << sumYRas << "\n";
    std::cout << "k: " << k << "\n";
    std::cout << "b: " << b << "\n";
    std::cout << "a0: " << a0 << "\n";
    std::cout << "a1: " << a1 << "\n";
    std::cout << "a2: " << a2 << "\n";

    return 0;
}

*/
