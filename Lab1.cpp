#include <iostream>
#include <fstream>
#include<cmath>
#include <iomanip>
#include<string>
#include<vector>
double k = 2.0;
double m = 0.3;
int n = 2;
int indp = 1;
double tau = 0.01;
double Tmax = 50.0;
double delta = 0.000001;
double eps = tau / 10;
int indikG = 0;
int indcorStep = 0;
int sizeB = 4;
double eps1 = 0.0000001;
std::vector<double> setka(double L1, double L2, int n)
{
    std::vector<double> points;
    double h1, h2;
    h1 = 2 * L1 / n;
    h2 = 2 * L2 / n;
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < n + 1; j++)
        {
            points.push_back(-L1 + i * h1);
            points.push_back(-L2 + j * h2);
        }
    }
    return points;
}
double AD(double (*funcs)(double, double*, int), double t, double* u, double delta, int i)
{
    double* x;
    x = new double[n];
    for (int k = 0; k < n; k++)
    {
        x[k] = u[k];
    }
    x[i] = x[i] + delta;
    return (funcs(t, x, n) - funcs(t, u, n)) / delta;
    delete[] x;
}
void SumVectCop(double* x, double* y, int n)
{
    for (int i = 0; i < n; i++)
    {
        x[i] = x[i] + y[i];
    }
}
void SumVect(double* x,double* y,double* result, int n)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = x[i] + y[i];
    }
}
void SclVect(double x, double* y, double* result, int n)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = x*y[i];
    }
}
void SclVectCop(double x, double* y, int n)
{
    for (int i = 0; i < n; i++)
    {
        y[i] = x * y[i];
    }
}
double NormVect(double* x, int n)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
        res = res + x[i] * x[i];
    }
    res = sqrt(res);
    return res;
}
double NormVect2(double* x,double* y, int n)
{
    double res = 0;
    for (int i = 0; i < n; i++)
    {
        res = res + (x[i]-y[i]) * (x[i]-y[i]);
    }
    res = sqrt(res);
    return res;
}
void MultWV(double** Matrix, double* x, double* result, int n)
{
    for (int i = 0; i < n; i++)
    {
        result[i] = 0;
        for (int j = 0; j < n; j++)
        {
            result[i] = result[i] + Matrix[i][j] * x[j];
        }
    }
}
int FindLeader(double** Matrix, int n, int j)
{
    int maxint = j;
    double max = 0;
    for (int i = j; i < n; i++)
    {
        if (abs(Matrix[i][j]) > max)
        {
            max = abs(Matrix[i][j]);
            maxint = i;
        }
    }
    if (max == 0)
    {
        indikG = 1;
    }
    return maxint;
}
void SwichLines(double** Matrix, double* rightb, int i, int j)
{
    double resb = rightb[i];
    double* resline = Matrix[i];
    Matrix[i] = Matrix[j];
    Matrix[j] = resline;
    rightb[i] = rightb[j];
    rightb[j] = resb;

}
void Gauss(double** Matrix, double* rightb, int n, double* x)
{
    double d, s;
    for (int k = 0; k < n - 1; k++)
    {
        SwichLines(Matrix, rightb, FindLeader(Matrix, n, k), k);
        if (indikG == 1)
        {
            return;
        }
        for (int j = k + 1; j < n; j++)
        {
            if (abs(Matrix[j][k]) != 0)
            {
                d = Matrix[j][k] / Matrix[k][k];
                Matrix[j][k] = 0;
                for (int i = k + 1; i < n; i++)
                {
                    Matrix[j][i] = Matrix[j][i] - d * Matrix[k][i];
                }
                rightb[j] = rightb[j] - d * rightb[k];
            }
        }
    }
    if (abs(Matrix[n - 1][n - 1]) < 0.00000000000001)
    {
        indikG = 1;
        return;
    }
    for (int k = n; k > 0; k--)
    {
        d = 0;
        for (int j = k; j < n; j++)
        {
            s = Matrix[k - 1][j] * x[j];
            d = d + s;
        }
        x[k - 1] = (rightb[k - 1] - d) / Matrix[k - 1][k - 1];
    }
}
double Func1(double t, double* u, int n)
{
    if (indp == 1)
    {
        return u[1];
    }
    if (indp == 2)
    {
        return u[1];
    }
    if (indp == 3)
    {
        return 0.4 * u[0] - 0.05 * u[0] * u[0] - 0.1 * u[0] * u[1];
    }
    
}
double Func2(double t, double* u, int n)
{
    if (indp == 1)
    {
        return -(k / m) * u[0];
    }
    if (indp == 2)
    {
        return 0.16 * cos(0.833 * t) - 0.05 * u[1] + 0.5 * u[0] - 0.5 * u[0] * u[0] * u[0];
    }
    if (indp == 3)
    {
        return -0.1 * u[1] - 0.003 * u[1] * u[1] + 0.08 * u[0] * u[1];
    }
}
double (*funcs[])(double , double* , int ) = { Func1,Func2 };
void newtonNEil(double (*funcs[])(double, double*, int), double t, double* u, double eps,double step, double* result)
{
    double** Matrix;
    double* rightb;
    double* res1;
    double* res2;
    double* dx;
    double* x;
    double* y;
    double tt = t + step;
    rightb = new double[n];
    dx = new double[n];
    x = new double[n];
    y = new double[n];
    res1 = new double[n];
    res2 = new double[n];
    Matrix = new double*[n];
    for (int i = 0; i < n; i++)
    {
        res1[i] = u[i];
        res2[i] = u[i] + 1.0;
        Matrix[i] = new double[n];
        rightb[i] = -step * funcs[i](tt, res1, n);
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                Matrix[i][j] = step * AD(funcs[i], tt, res1, delta, j) - 1.0;
            }
            else
            {
                Matrix[i][j] = step * AD(funcs[i], tt, res1, delta, j);
            }
        }
    }
    Gauss(Matrix, rightb, n, dx);
    SumVect(res1, dx, res2, n);
    while (NormVect2(res1,res2,n)>eps)
    {
        
        for (int i = 0; i < n; i++)
        {
            res1[i] = res2[i];
        }
        for (int i = 0; i < n; i++)
        {
            x[i] = -u[i];
            dx[i] = -step * funcs[i](tt, res1, n);
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    Matrix[i][j] = step * AD(funcs[i], tt, res1, delta, j) - 1.0;
                }
                else
                {
                    Matrix[i][j] = step * AD(funcs[i], tt, res1, delta, j);
                }
            }
        }
        SumVect(dx, x, y, n);
        SumVect(y, res1, rightb, n);
        Gauss(Matrix, rightb, n, dx);
        SumVect(res1, dx, res2, n);
    }
    for (int i = 0; i < n; i++)
    {
        result[i] = res2[i];
        delete[] Matrix[i];
    }
    delete[] rightb;
    delete[] dx;
    delete[] x;
    delete[] y;
    delete[] res1;
    delete[] res2;
    delete[] Matrix;
}
void newtonSimetr(double (*funcs[])(double, double*, int), double t, double* u, double eps, double step, double* result)
{
    double** Matrix;
    double* rightb;
    double* res1;
    double* res2;
    double* dx;
    double* x;
    double* y;
    double tt = t + step;
    rightb = new double[n];
    dx = new double[n];
    x = new double[n];
    y = new double[n];
    res1 = new double[n];
    res2 = new double[n];
    Matrix = new double* [n];
    for (int i = 0; i < n; i++)
    {
        res1[i] = u[i];
        res2[i] = u[i] + 1.0;
        Matrix[i] = new double[n];
        rightb[i] = -step/2 * (funcs[i](tt, res1, n)+ funcs[i](t, u, n));
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                Matrix[i][j] = step /2 * AD(funcs[i], tt, res1, delta, j) - 1.0;
            }
            else
            {
                Matrix[i][j] = step /2 * AD(funcs[i], tt, res1, delta, j);
            }
        }
    }
    Gauss(Matrix, rightb, n, dx);
    SumVect(res1, dx, res2, n);
    while (NormVect2(res1, res2, n) > eps)
    {
        for (int i = 0; i < n; i++)
        {
            res1[i] = res2[i];
        }
        for (int i = 0; i < n; i++)
        {
            x[i] = -u[i];
            dx[i] = -step /2 * (funcs[i](tt, res2, n)+ funcs[i](t, u, n));
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    Matrix[i][j] = step /2 * AD(funcs[i], tt, res2, delta, j) - 1.0;
                }
                else
                {
                    Matrix[i][j] = step /2 * AD(funcs[i], tt, res2, delta, j);
                }
            }
        }
        SumVect(dx, x, y, n);
        SumVect(y, res1, rightb, n);
        Gauss(Matrix, rightb, n, dx);
        SumVect(res1, dx, res2, n);
    }
    for (int i = 0; i < n; i++)
    {
        result[i] = res2[i];
        delete[] Matrix[i];
    }
    delete[] rightb;
    delete[] dx;
    delete[] x;
    delete[] y;
    delete[] res1;
    delete[] res2;
    delete[] Matrix;
}
void YavnEilerStep(double** res, double t,double step, int n)
{
    double* x;
    double* y;
    double* z;
    x = new double[n];
    y = new double[n];
    z = new double[n];
    for (int i = 0; i < n; i++)
    {
        y[i] = res[0][i];
        x[i] = funcs[i](t, res[0], n);
    }
    SclVect(step, x, z, n);
    SumVect(y, z, res[1], n);
    delete[] x;
    delete[] y;
    delete[] z;
}
void YavnEiler(double t0, double* u0, int n, std::string nn)
{
    double tautek = tau;
    double** res;
    res = new double* [2];
    for (int i = 0; i < 2; i++)
    {
        res[i] = new double[n];
    }
    double** res1;
    double contr;
    res1 = new double* [2];
    for (int i = 0; i < 2; i++)
    {
        res1[i] = new double[n];
    }
    double t = t0;
    for (int i = 0; i < n; i++)
    {
        res[0][i] = u0[i];
    }
    std::ofstream ans;
    std::string name;
    name = "YEil";
    name = nn + name;
    if (indp == 1)
    {
        name = name + "test1";
    }
    if (indp == 2)
    {
        name = name + "v16";
    }
    if (indp == 3)
    {
        name = name + "v20";
    }
    name = name + ".txt";
    ans.open(name);
    ans << std::setprecision(15);
    ans << t << ' ';
    for (int i = 0; i < n; i++)
    {
        ans << res[0][i] << ' ';
    }
    ans << std::endl;
    while (t <= Tmax)
    {
        
        if (indcorStep == 0)
        {
            YavnEilerStep(res, t, tautek, n);
            t = t + tautek;
            for (int i = 0; i < n; i++)
            {
                res[0][i] = res[1][i];
            }
            ans << t << ' ';
            for (int i = 0; i < n; i++)
            {
                ans << res[0][i] << ' ';
            }
            ans << std::endl;
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                
                res1[0][i] = res[0][i];
            }
            YavnEilerStep(res, t, tautek, n);
            YavnEilerStep(res1, t, tautek/2, n);
            for (int i = 0; i < n; i++)
            {
                res1[0][i] = res1[1][i];
            }
            YavnEilerStep(res1, t+tautek/2, tautek/2, n);
            contr= ((1. / 3.) * NormVect2(res[1], res1[1], n));
            if (contr > eps1)
            {
                tautek = tautek / 2;
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    res[0][i] = res1[1][i];
                }
                t = t + tautek;
                ans << t << ' ';
                for (int i = 0; i < n; i++)
                {
                    ans << res[0][i] << ' ';
                }
                ans << std::endl;
                if (contr < (10 * eps1))
                {
                    tautek = tautek * 2;
                }
            }
        }
        
    }
    ans.close();
    for (int i = 0; i < 2; i++)
    {
        delete[] res[i];
        delete[] res1[i];
    }
    delete[] res;
    delete[] res1;
}
void NeyavnEilerStep(double** res, double t, double step, int n)
{
    eps = tau / 100;
    newtonNEil(funcs, t, res[0], eps, step, res[1]);
}
void NeyavnEiler(double t0, double* u0, int n, std::string nn)
{
    double** res;
    res = new double* [2];
    for (int i = 0; i < 2; i++)
    {
        res[i] = new double[n];
    }
    double tautek = tau;
    double** res1;
    double contr;
    res1 = new double* [2];
    for (int i = 0; i < 2; i++)
    {
        res1[i] = new double[n];
    }
    double t = t0;
    for (int i = 0; i < n; i++)
    {
        res[0][i] = u0[i];
    }
    std::ofstream ans;
    std::string name;
    name = "NeyEil";
    name = nn + name;
    if (indp == 1)
    {
        name = name + "test1";
    }
    if (indp == 2)
    {
        name = name + "v16";
    }
    if (indp == 3)
    {
        name = name + "v20";
    }
    name = name + ".txt";
    ans.open(name);
    ans << std::setprecision(15);
    ans << t << ' ';
    for (int i = 0; i < n; i++)
    {
        ans << res[0][i] << ' ';
    }
    ans << std::endl;
    while (t <= Tmax)
    {

        if (indcorStep == 0)
        {
            NeyavnEilerStep(res, t, tautek, n);
            t = t + tautek;
            for (int i = 0; i < n; i++)
            {
                res[0][i] = res[1][i];
            }
            ans << t << ' ';
            for (int i = 0; i < n; i++)
            {
                ans << res[0][i] << ' ';
            }
            ans << std::endl;
        }
        else
        {
            for (int i = 0; i < n; i++)
            {

                res1[0][i] = res[0][i];
            }
            NeyavnEilerStep(res, t, tautek, n);
            NeyavnEilerStep(res1, t, tautek / 2, n);
            for (int i = 0; i < n; i++)
            {
                res1[0][i] = res1[1][i];
            }
            NeyavnEilerStep(res1, t + tautek / 2, tautek / 2, n);
            contr = ((1. / 3.) * NormVect2(res[1], res1[1], n));
            if (contr > eps1)
            {
                tautek = tautek / 2;
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    res[0][i] = res1[1][i];
                }
                t = t + tautek;
                ans << t << ' ';
                for (int i = 0; i < n; i++)
                {
                    ans << res[0][i] << ' ';
                }
                ans << std::endl;
                if (contr < (10 * eps1))
                {
                    tautek = tautek * 2;
                }
            }
        }

    }
    ans.close();
    for (int i = 0; i < 2; i++)
    {
        delete[] res[i];
        delete[] res1[i];
    }
    delete[] res;
    delete[] res1;
}
void SimetrStep(double** res, double t, double step, int n)
{
    eps = tau / 100;
    newtonSimetr(funcs, t, res[0], eps, step, res[1]);
}
void Simetr(double t0, double* u0, int n, std::string nn)
{
    double** res;
    res = new double* [2];
    for (int i = 0; i < 2; i++)
    {
        res[i] = new double[n];
    }
    double tautek = tau;
    double** res1;
    double contr;
    res1 = new double* [2];
    for (int i = 0; i < 2; i++)
    {
        res1[i] = new double[n];
    }
    double t = t0;
    for (int i = 0; i < n; i++)
    {
        res[0][i] = u0[i];
    }
    std::ofstream ans;
    std::string name;
    name = "Simetr";
    name = nn + name;
    if (indp == 1)
    {
        name = name + "test1";
    }
    if (indp == 2)
    {
        name = name + "v16";
    }
    if (indp == 3)
    {
        name = name + "v20";
    }
    name = name + ".txt";
    ans.open(name);
    ans << std::setprecision(15);
    ans << t << ' ';
    for (int i = 0; i < n; i++)
    {
        ans << res[0][i] << ' ';
    }
    ans << std::endl;
    while (t <= Tmax)
    {

        if (indcorStep == 0)
        {
            SimetrStep(res, t, tautek, n);
            t = t + tautek;
            for (int i = 0; i < n; i++)
            {
                res[0][i] = res[1][i];
            }
            ans << t << ' ';
            for (int i = 0; i < n; i++)
            {
                ans << res[0][i] << ' ';
            }
            ans << std::endl;
        }
        else
        {
            for (int i = 0; i < n; i++)
            {

                res1[0][i] = res[0][i];
            }
            SimetrStep(res, t, tautek, n);
            SimetrStep(res1, t, tautek / 2, n);
            for (int i = 0; i < n; i++)
            {
                res1[0][i] = res1[1][i];
            }
            SimetrStep(res1, t + tautek / 2, tautek / 2, n);
            contr = ((1. / 3.) * NormVect2(res[1], res1[1], n));
            if (contr > eps1)
            {
                tautek = tautek / 2;
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    res[0][i] = res1[1][i];
                }
                t = t + tautek;
                ans << t << ' ';
                for (int i = 0; i < n; i++)
                {
                    ans << res[0][i] << ' ';
                }
                ans << std::endl;
                if (contr < (10 * eps1))
                {
                    tautek = tautek * 2;
                }
            }
        }

    }
    ans.close();
    for (int i = 0; i < 2; i++)
    {
        delete[] res[i];
        delete[] res1[i];
    }
    delete[] res;
    delete[] res1;
    
}
void RungeKuttStep(double** res, double t, double step, int n, double** Butch)
{
    double tt;
    double* y;
    double* z;
    double** kk;
    kk = new double* [sizeB];
    y = new double[n];
    z = new double[n];
    for (int i = 0; i < sizeB; i++)
    {
        kk[i] = new double[n];
    }
    for (int i = 0; i < n; i++)
    {
        res[1][i] = res[0][i];
    }
    for (int i = 0; i < n; i++)
    {
        kk[0][i] = funcs[i](t, res[1], n);
    }
    for (int i = 0; i < sizeB - 1; i++)
    {
        for (int p = 0; p < n; p++)
        {
            y[p] = 0.0;
        }
        for (int j = 1; j < sizeB; j++)
        {
            SclVect(Butch[i][j], kk[j-1], z, n);
            SumVectCop(y, z, n);
        }
        SclVectCop(step, y, n);
        SumVectCop(y, res[1], n);
        tt = t + Butch[i][0] * step;
        for (int m = 0; m < n; m++)
        {
            kk[i + 1][m] = funcs[m](tt, y, n);
        }
    }
    for (int i = 0; i < sizeB; i++)
    {
        SclVectCop(Butch[sizeB - 1][i], kk[i], n);
        SclVectCop(step, kk[i], n);
        SumVectCop(res[1], kk[i], n);
    }
    for (int i = 0; i < sizeB; i++)
    {
        delete[] kk[i];
    }
    delete[] y;
    delete[] z;
    delete[] kk;
}
void RungeKutt(double t0, double* u0, int n, double** Butch, std::string nn)
{
    double** res;
    double** res1;
    double contr;
    res1 = new double* [2];
    res = new double* [2];
    for (int i = 0; i < 2; i++)
    {
        res[i] = new double[n];
        res1[i] = new double[n];
    }
    double tautek = tau;
    double t = t0;
    for (int i = 0; i < n; i++)
    {
        res[0][i] = u0[i];
    }
    std::ofstream ans;
    std::string name;
    name = "RungKutt";
    name = nn + name;
    if (indp == 1)
    {
        name = name + "test1";
    }
    if (indp == 2)
    {
        name = name + "v16";
    }
    if (indp == 3)
    {
        name = name + "v20";
    }
    name = name + ".txt";
    ans.open(name);
    ans << std::setprecision(15);
    ans << t << ' ';
    for (int i = 0; i < n; i++)
    {
        ans << res[0][i] << ' ';
    }
    ans << tautek << ' ' << std::cos(std::sqrt(k / m) * t) << std::endl;
    while (t <= Tmax)
    {

        if (indcorStep == 0)
        {
            RungeKuttStep(res, t, tautek, n, Butch);
            t = t + tautek;
            for (int i = 0; i < n; i++)
            {
                res[0][i] = res[1][i];
            }
            ans << t << ' ';
            for (int i = 0; i < n; i++)
            {
                ans << res[0][i] << ' ';
            }
            ans << std::endl;
        }
        else
        {
            for (int i = 0; i < n; i++)
            {

                res1[0][i] = res[0][i];
            }
            RungeKuttStep(res, t, tautek, n, Butch);
            RungeKuttStep(res1, t, tautek / 2, n, Butch);
            for (int i = 0; i < n; i++)
            {
                res1[0][i] = res1[1][i];
            }
            RungeKuttStep(res1, t + tautek / 2, tautek / 2, n, Butch);
            contr = ((1. / 15.) * NormVect2(res[1], res1[1], n));
            if (contr > eps1)
            {
                tautek = tautek / 2;
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    res[0][i] = res1[1][i];
                }
                t = t + tautek;
                ans << t << ' ';
                for (int i = 0; i < n; i++)
                {
                    ans << res[0][i] << ' ';
                }
                ans << tautek << ' ' << std::cos(std::sqrt(k / m) * t) << std::endl;
                if (contr < (10 * eps1))
                {
                    tautek = tautek * 2;
                }
            }
        }

    }
    ans.close();
    for (int i = 0; i < 2; i++)
    {
        delete[] res[i];
        delete[] res1[i];
    }
    delete[] res;
    delete[] res1;
    

}
void AdamsBoshStep(double** res,double** f, double t, double step, int n, double* CB)
{
    double* y;
    y = new double[n];
    double* resf;
    resf = new double[n];
    for (int i = 0; i < n; i++)
    {
        res[1][i]=res[0][i];
    }
    for (int i = 0; i < n; i++)
    {
        resf[i] = funcs[i](t, res[1], n);
    }
    for (int i = 1; i < 4; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f[i - 1][j] = f[i][j];
        }
    }
    for (int j = 0; j < n; j++)
    {
        f[3][j] = resf[j];
    }
    for (int i = 0; i < 4; i++)
    {

        SclVect(CB[i], f[i], y, n);
        SclVectCop(step, y, n);
        SumVectCop(res[1], y, n);
    }
    delete[] y;
    delete[] resf;
}
void AdamsBoshf(double t0, double* u0, int n, double** Butch, double* CB,std::string nn)
{
        double t = t0;
        double** f;
        f = new double*[4];
        double** res;
        res = new double* [2];
        for (int i = 0; i < 2; i++)
        {
            res[i] = new double[n];
        }
        for (int i = 0; i < 4; i++)
        {
            f[i] = new double[n];
        }
        for (int i = 0; i < n; i++)
        {
            res[0][i] = u0[i];
        }
        std::ofstream ans;
        std::string name;
        name = "AdamsBoshf";
        name = nn+name;
        if (indp == 1)
        {
            name = name + "test1";
        }
        if (indp == 2)
        {
            name = name + "v16";
        }
        if (indp == 3)
        {
            name = name + "v20";
        }
        name = name + ".txt";
        ans.open(name);
        ans << std::setprecision(15);
        for(int q = 0; q < 4;q++)
        {
            for (int i = 0; i < n; i++)
            {
                f[q][i]= funcs[i](t, res[0], n);
            }
            ans << t << ' ';
            for (int i = 0; i < n; i++)
            {
                ans << res[0][i] << ' ';
            }
            ans << std::endl;
            RungeKuttStep(res, t, tau, n, Butch);
            for (int i = 0; i < n; i++)
            {
                res[0][i] = res[1][i];
            }
            t = t + tau;
        }
        while (t <= Tmax)
        {
            ans << t << ' ';
            for (int i = 0; i < n; i++)
            {
                ans << res[0][i] << ' ';
            }
            ans << std::endl;
            AdamsBoshStep(res, f, t, tau, n, CB);
            for (int i = 0; i < n; i++)
            {
                res[0][i] = res[1][i];
            }
            t = t + tau;
        }
        ans.close();
        for (int i = 0; i < 4; i++)
        {
            delete[] f[i];
        }
        for (int i = 0; i < 2; i++)
        {
            delete[] res[i];
        }
        delete[] res;
        delete[] f;
}
void PrognCorStep(double** res, double** f, double t, double step, int n, double* CB, double* CBcor)
{
    double* y;
    y = new double[n];
    double* x;
    x = new double[n];
    double* resf;
    double** f0;
    f0 = new double* [4];
    resf = new double[n];
    for (int i = 0; i < n; i++)
    {
        res[1][i] = res[0][i];
    }
    for (int i = 0; i < 4; i++)
    {
        f0[i] = new double[n];
    }
    for (int i = 0; i < n; i++)
    {
        resf[i] = funcs[i](t, res[0], n);
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = res[0][i];
    }
    for (int i = 1; i < 4; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f[i - 1][j] = f[i][j];
        }
    }
    for (int j = 0; j < n; j++)
    {
        f[3][j] = resf[j];
    }
    for (int i = 0; i < 4; i++)
    {

        SclVect(CB[i], f[i], y, n);
        SclVectCop(step, y, n);
        SumVectCop(x, y, n);
    }
    for (int i = 0; i < n; i++)
    {
        resf[i] = funcs[i](t + step, x, n);
    }
    for (int i = 1; i < 4; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f0[i - 1][j] = f[i][j];
        }
    }
    for (int j = 0; j < n; j++)
    {
        f0[3][j] = resf[j];
    }
    for (int i = 0; i < 4; i++)
    {

        SclVect(CBcor[i], f0[i], y, n);
        SclVectCop(step, y, n);
        SumVectCop(res[1], y, n);
    }
    for (int i = 0; i < 4; i++)
    {
        delete[] f0[i];
    }
    delete[] f0;
    delete[] x;
    delete[] y;
    delete[] resf;
}
void PrognCor(double t0, double* u0, int n, double** Butch, double* CB,double* CBcor, std::string nn)
{
    double t = t0;
    double** f;
    f = new double* [4];
    double** res;
    res = new double* [2];
    for (int i = 0; i < 2; i++)
    {
        res[i] = new double[n];
    }
    for (int i = 0; i < 4; i++)
    {
        f[i] = new double[n];
    }
    for (int i = 0; i < n; i++)
    {
        res[0][i] = u0[i];
    }
    std::ofstream ans;
    std::string name;
    name = "PrognCor";
    name = nn + name;
    if (indp == 1)
    {
        name = name + "test1";
    }
    if (indp == 2)
    {
        name = name + "v16";
    }
    if (indp == 3)
    {
        name = name + "v20";
    }
    name = name + ".txt";
    ans.open(name);
    ans << std::setprecision(15);
    for (int q = 0; q < 4; q++)
    {
        for (int i = 0; i < n; i++)
        {
            f[q][i] = funcs[i](t, res[0], n);
        }
        ans << t << ' ';
        for (int i = 0; i < n; i++)
        {
            ans << res[0][i] << ' ';
        }
        ans << std::endl;
        RungeKuttStep(res, t, tau, n, Butch);
        for (int i = 0; i < n; i++)
        {
            res[0][i] = res[1][i];
        }
        t = t + tau;
    }
    while (t <= Tmax)
    {
        for (int i = 0; i < n; i++)
        {
            res[0][i] = res[1][i];
        }
        ans << t << ' ';
        for (int i = 0; i < n; i++)
        {
            ans << res[0][i] << ' ';
        }
        ans << std::endl;
        PrognCorStep(res, f, t, tau, n, CB, CBcor);
        t = t + tau;
    }
    ans.close();
    for (int i = 0; i < 4; i++)
    {
        delete[] f[i];
    }
    for (int i = 0; i < 2; i++)
    {
        delete[] res[i];
    }
    delete[] res;
    delete[] f;
}
void test2()
{
    double** res;
    res = new double* [2];
    for (int i = 0; i < 2; i++)
    {
        res[i] = new double[n];
    }
    double L1, L2;
    int m=50;
    L1 = 2;
    L2 = 2;
    std::string nn="";
    std::vector<double> p;
    p = setka(L1, L2, m);
    std::ofstream ans;
    ans.open("Phaz2.txt");
    double** Butch;
    double* CB;
    double* CBcor;
    double t0=0;
    double* u0;
    double t = t0;
    u0 = new double[n];
    CB = new double[4];
    CBcor = new double[4];
    Butch = new double* [sizeB];
    for (int i = 0; i < sizeB; i++)
    {
        Butch[i] = new double[sizeB];
        for (int j = 0; j < sizeB; j++)
        {
            Butch[i][j] = 0.0;
        }
    }
    CBcor[0] = (1.0) / (24.0);
    CBcor[1] = -(5.0) / (24.0);
    CBcor[2] = (19.0) / (24.0);
    CBcor[3] = (9.0) / (24.0);
    CB[0] = -(9.0) / (24.0);
    CB[1] = (37.0) / (24.0);
    CB[2] = -(59.0) / (24.0);
    CB[3] = (55.0) / (24.0);
    Butch[0][0] = (1.0) / (2.0);
    Butch[0][1] = (1.0) / (2.0);
    Butch[1][0] = (1.0) / (2.0);
    Butch[1][2] = (1.0) / (2.0);
    Butch[2][0] = 1.0;
    Butch[2][3] = 1.0;
    Butch[3][0] = (1.0) / (6.0);
    Butch[3][1] = (1.0) / (3.0);
    Butch[3][2] = (1.0) / (3.0);
    Butch[3][3] = (1.0) / (6.0);
    nn = "";
    tau = 0.1;
    indp = 2;
    indcorStep = 0;
    eps1 = 0.000001;
    if (indp == 1)
    {
        t0 = 0;
        Tmax = 5.;
        u0[0] = 0.;
        u0[1] = 1.;
    }
    if (indp == 2)
    {
        t0 = 0;
        Tmax = 50.0;
        u0[0] = 0.1;
        u0[1] = -0.1;
    }
    if (indp == 3)
    {
        t0 = 0;
        Tmax = 150.0;
        u0[0] = 1.;
        u0[1] = 4.;
    }
    t = t0;
    while (t < 8)
    {
        for (int i = 0; i < p.size(); i++)
        {
            res[0][0] = p[i];
            res[0][1] = p[i + 1];
            RungeKuttStep(res, t, tau, n, Butch);
            ans << t << ' ' << res[0][0] << ' ' << res[0][1] << ' ' << res[1][0] << ' ' << res[1][1] << std::endl;
            i++;
        }
        t = t + 0.05;
    }
    ans.close();
    for (int i = 0; i < 2; i++)
    {
        delete[] res[i];
    }
    delete[] res;
    
}
int main()
{
    std::string nn;
    double** Butch;
    double* CB;
    double* CBcor;
    double t0;
    double* u0;
    u0 = new double[n];
    CB = new double[4];
    CBcor = new double[4];
    Butch = new double* [sizeB];
    for (int i = 0; i < sizeB; i++)
    {
        Butch[i] = new double[sizeB];
        for (int j = 0; j < sizeB; j++)
        {
            Butch[i][j] = 0.0;
        }
    }
    CBcor[0] = (1.0) / (24.0);
    CBcor[1] = -(5.0) / (24.0);
    CBcor[2] = (19.0) / (24.0);
    CBcor[3] = (9.0) / (24.0);
    CB[0] = -(9.0) / (24.0);
    CB[1] = (37.0) / (24.0);
    CB[2] = -(59.0) / (24.0);
    CB[3] = (55.0) / (24.0);
    Butch[0][0] = (1.0) / (2.0);
    Butch[0][1] = (1.0) / (2.0);
    Butch[1][0] = (1.0) / (2.0);
    Butch[1][2] = (1.0) / (2.0);
    Butch[2][0] = 1.0;
    Butch[2][3] = 1.0;
    Butch[3][0] = (1.0) / (6.0);
    Butch[3][1] = (1.0) / (3.0);
    Butch[3][2] = (1.0) / (3.0);
    Butch[3][3] = (1.0) / (6.0);
    nn = "0.00078125";
    tau = 0.00078125;
    indp = 2;
    indcorStep = 0;
    eps1 = 0.000001;
    if (indp == 1)
    {
        t0 = 0;
        Tmax = 50.;
        u0[0] = 0.;
        u0[1] = 1.;
    }
    if (indp == 2)
    {
        t0 = 0;
        Tmax = 100.0;
        u0[0] = 0.1;
        u0[1] = -0.1;
    }
    if (indp == 3)
    {
        t0 = 0;
        Tmax = 150.0;
        u0[0] = 1.;
        u0[1] = 4.;
    }
   // test2();
    //nn = "fgr";
    //YavnEiler(t0, u0, n,nn);
    //NeyavnEiler(t0, u0, n,nn);
   // Simetr(t0, u0, n,nn);
    RungeKutt(t0, u0, n, Butch,nn);
    AdamsBoshf(t0, u0, n, Butch, CB,nn);
    PrognCor(t0, u0, n, Butch, CB, CBcor,nn);
    std::cout << "Hello World!\n";
    
    delete[] u0;
    
}

