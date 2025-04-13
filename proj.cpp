#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <string>

#define N_ABCD 951
#define A1B 7
#define A2 -1
#define A3 -1
#define F 4
#define MIN_ERR 1.0e-9
#define A1C 3
#define ERR_TAB_SIZE 3000
#define TIME_FILE_NAME "TimeE.txt"

using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;
using std::chrono::milliseconds;
using std::chrono::time_point;

// tu są zaimplemenowane metody podane w poleceniu

// stopnie = radiany * (180 / pi)
// radiany = stopnie / (180/pi)

void saveErrToFile(double *v, string fileName, int N)
{
    ofstream file;
    file.open(fileName, std::ios_base::openmode::_S_trunc);
    for (int i = 0; i < N; i++)
    {
        file << v[i] << endl;
    }
    file.close();
}

void saveTimeToFile(int time)
{
    ofstream file;
    file.open(TIME_FILE_NAME, std::ios::app);
    file << time << endl;
    file.close();
}

void setUpDataBE(double **A, double *b, double *x, int N)
{
    for (int i = 0; i < N; i++)
    {
        b[i] = sin((i * F) / (180 / M_PI));
        x[i] = 1;
        for (int j = 0; j < N; j++)
        {
            if (i == j)
                A[i][j] = A1B;
            else if (i == j + 1 || j == i + 1)
                A[i][j] = A2;
            else if (j == i + 2 || i == j + 2)
                A[i][j] = A3;
            else
                A[i][j] = 0;
        }
    }
}

void setUpDataCD(double **A, double *b, double *x, int N)
{
    for (int i = 0; i < N; i++)
    {
        b[i] = sin((i * F) / (180 / M_PI));
        x[i] = 1;
        for (int j = 0; j < N; j++)
        {
            if (i == j)
                A[i][j] = A1C;
            else if (i == j + 1 || j == i + 1)
                A[i][j] = A2;
            else if (j == i + 2 || i == j + 2)
                A[i][j] = A3;
            else
                A[i][j] = 0;
        }
    }
}

double eucNorm(double *v, int N)
{
    double norm = 0;

    for (int i = 0; i < N; i++)
        norm += pow(v[i], 2);

    return sqrt(norm);
}

void mulitplyByVector(double **A, double *v, double *result, int N)
{
    for (int i = 0; i < N; i++)
    {
        result[i] = 0.0;
        for (int j = 0; j < N; j++)
            result[i] += A[i][j] * v[j];
    }
}

void Jacobi(double **A, double *b, double *x, int N, bool save = 0, string fileName = "", bool saveTime = 0)
{
    int iterations = 0;
    double norm;
    double *r = new double[N];
    double *err = new double[ERR_TAB_SIZE];
    // x - prev
    // r - newResult

    time_point start = high_resolution_clock::now();

    do
    {
        double element = 0.0;
        for (int i = 0; i < N; i++)
        {
            element = 0.0;
            for (int j = 0; j < N; j++)
                if (j != i)
                    element += A[i][j] * x[j];

            r[i] = (b[i] - element) / A[i][i]; // r obecny wynik
        }

        for (int i = 0; i < N; i++)
            x[i] = r[i]; // obecny wynik ustlany jako poprzedni wynik w nastepnej iteracji
        // r == x

        mulitplyByVector(A, x, r, N);
        for (int i = 0; i < N; i++)
            r[i] -= b[i];
        norm = eucNorm(r, N);
        err[iterations] = norm;
        iterations++;
    } while (norm > MIN_ERR);

    time_point stop = high_resolution_clock::now();
    microseconds duration = duration_cast<milliseconds>(stop - start);
    cout << "Jacobi method" << endl;
    cout << "Time in miliseconds: " << duration.count() / 1000 << endl;
    cout << "Iterations: " << iterations << endl;
    cout << "Norm: " << norm << endl;

    if (save)
        saveErrToFile(err, fileName, iterations);

    if (saveTime)
        saveTimeToFile(duration.count() / 1000);

    delete r;
    delete err;
}

void GaussSeidel(double **A, double *b, double *x, int N, bool save = 0, string fileName = "", bool saveTime = 0)
{
    int iterations = 0;
    double norm;
    double *r = new double[N];
    double *err = new double[ERR_TAB_SIZE];
    for (int i = 0; i < N; i++)
        r[i] = 1.0;
    // x - prev (k)
    // r - newResult (k + 1)

    time_point start = high_resolution_clock::now();

    do
    {
        double element = 0.0;
        for (int i = 0; i < N; i++)
        {
            element = 0.0;
            for (int j = 0; j < i; j++) // górna część (k+1)
                element += A[i][j] * r[j];
            for (int j = i + 1; j < N; j++) // dolna część (k)
                element += A[i][j] * x[j];

            r[i] = (b[i] - element) / A[i][i]; // r obecny wynik
        }

        for (int i = 0; i < N; i++)
            x[i] = r[i]; // obecny wynik ustlany jako poprzedni wynik w nastepnej iteracji
        // r == x

        mulitplyByVector(A, x, r, N);
        for (int i = 0; i < N; i++) // do obliczenia błędu
            r[i] -= b[i];

        norm = eucNorm(r, N);
        err[iterations] = norm;
        iterations++;

        for (int i = 0; i < N; i++) // odwracanie aby wyniki były poprawne
        {
            r[i] += b[i];
        }
    } while (norm > MIN_ERR);

    time_point stop = high_resolution_clock::now();
    microseconds duration = duration_cast<milliseconds>(stop - start);
    cout << "Gauss-Seidel method" << endl;
    cout << "Time in miliseconds: " << duration.count() / 1000 << endl;
    cout << "Iterations: " << iterations << endl;
    cout << "Norm: " << norm << endl;

    if (save)
        saveErrToFile(err, fileName, iterations);
    if (saveTime)
        saveTimeToFile(duration.count() / 1000);

    delete r;
    delete err;
}

void LU(double **A, double *b, double *x, int N, bool saveTime = 0)
{
    // setup
    double norm;
    double *r = new double[N];
    double **U = new double *[N]; // górna trójkątna część
    double **L = new double *[N]; // dolna trójkątna część
    double *y = new double[N];

    for (int i = 0; i < N; i++)
    {
        y[i] = 1.0;
        U[i] = new double[N];
        L[i] = new double[N];
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
                L[i][j] = 1;
            else
                L[i][j] = 0;
            U[i][j] = A[i][j];
        }
    }

    time_point start = high_resolution_clock::now();

    for (int i = 2; i < N + 1; i++)
    {
        for (int j = 1; j < i; j++)
        {
            L[i - 1][j - 1] = U[i - 1][j - 1] / U[j - 1][j - 1];
            for (int k = 0; k < N; k++) // k zastapuje numer kolumny
                U[i - 1][k] = U[i - 1][k] - L[i - 1][j - 1] * U[j - 1][k];
        }
    }

    for (int i = 0; i < N; i++) // podstawienie w przód
    {
        double element = 0.0;

        for (int j = 0; j < i; j++)
            element += L[i][j] * y[j];

        y[i] = (b[i] - element) / L[i][i];
    }

    for (int i = N - 1; i >= 0; i--) // podstawienie w tył
    {
        double element = 0.0;

        for (int j = i + 1; j < N; j++)
            element += U[i][j] * x[j];

        x[i] = (y[i] - element) / U[i][i];
    }
    mulitplyByVector(A, x, r, N);
    for (int i = 0; i < N; i++)
        r[i] -= b[i];

    norm = eucNorm(r, N);

    time_point stop = high_resolution_clock::now();
    microseconds duration = duration_cast<milliseconds>(stop - start);
    cout << "LU method" << endl;
    cout << "Time in miliseconds: " << duration.count() / 1000 << endl;
    cout << "Norm: " << norm << endl;

    if (saveTime)
        saveTimeToFile(duration.count() / 1000);

    delete r;
    delete y;
    for (int i = 0; i < N; i++)
    {
        delete U[i];
        delete L[i];
    }
    delete U;
    delete L;
}

int main()
{
    double **A = new double *[N_ABCD];
    for (int i = 0; i < N_ABCD; i++)
        A[i] = new double[N_ABCD];
    double *b = new double[N_ABCD];
    double *x = new double[N_ABCD];

    // wyczyść plik z czasem
    saveErrToFile(b,TIME_FILE_NAME,0);

    // B
    cout << endl;
    cout << "B: (N: " << N_ABCD << " A1: " << A1B << " A2,A3: " << A2 << " F: " << F << " )" << endl;
    setUpDataBE(A, b, x, N_ABCD);
    Jacobi(A, b, x, N_ABCD, 1, "JacobiErrB.txt");
    cout << endl;
    setUpDataBE(A, b, x, N_ABCD);
    GaussSeidel(A, b, x, N_ABCD, 1, "GaussSeidelB.txt");
    cout << endl;
    // C
    cout << endl;
    cout << "C: (N: " << N_ABCD << " A1: " << A1C << " A2,A3: " << A2 << " F: " << F << " )" << endl;
    setUpDataCD(A, b, x, N_ABCD);
    Jacobi(A, b, x, N_ABCD, 1, "JacobiErrC.txt");
    cout << endl;
    setUpDataCD(A, b, x, N_ABCD);
    GaussSeidel(A, b, x, N_ABCD, 1, "GaussSeidelC.txt");
    cout << endl;
    // D
    cout << endl;
    cout << "D: (N: " << N_ABCD << " A1: " << A1C << " A2,A3: " << A2 << " F: " << F << " )" << endl;
    setUpDataCD(A, b, x, N_ABCD);
    LU(A, b, x, N_ABCD);
    cout << endl;

    delete x;
    delete b;
    for (int i = 0; i < N_ABCD; i++)
        delete A[i];
    delete A;

    // E
    int sizes[] = {100, 500, 1000, 2000, 3000, 4000, 5000};
    cout << endl;
    cout << "E: (N: "
         << "100, 500, 1000, 2000, 3000, 4000, 5000"
         << " A1: " << A1B << " A2,A3: " << A2 << " F: " << F << " )" << endl;
    for (int i = 0; i < 7; i++)
    {
        double **A = new double *[sizes[i]];
        for (int j = 0; j < sizes[i]; j++)
            A[j] = new double[sizes[i]];
        double *b = new double[sizes[i]];
        double *x = new double[sizes[i]];
        cout << "Size: " << sizes[i] << endl;
        cout << endl;
        setUpDataBE(A, b, x, sizes[i]);
        Jacobi(A, b, x, sizes[i], 0, "", 1);
        cout << endl;
        setUpDataBE(A, b, x, sizes[i]);
        GaussSeidel(A, b, x, sizes[i], 0, "", 1);
        cout << endl;
        setUpDataBE(A, b, x, sizes[i]);
        LU(A, b, x, sizes[i], 1);
        cout << endl;

        delete x;
        delete b;
        for (int i = 0; i < sizes[i]; i++)
            delete A[i];
        delete A;
    }
    cout << endl;
}