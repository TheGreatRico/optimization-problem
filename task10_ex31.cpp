#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
using namespace std;

#define N 2       // Matrix dimensions
#define EPSILON 1e-10
const double tol = 1e-10;
const double beta = 10; // 0, 0.1, 1, 10
// initial system functions
double dot_x(double t, double x, double y, double px, double py) {
    return y;
}
double dot_y(double t, double x, double y, double px, double py) {
    return py * (1 + beta * x * x);
}
double dot_px(double t, double x, double y, double px, double py) {
    return beta * x * py * py;
}
double dot_py(double t, double x, double y, double px, double py) {
    return -px;
}
// R-K 5:
void shoot(double x0, double y0, double px0, double py0, double t0, double t_end, double ret_arr[]) { // h?
    double t_n = t0, x_n = x0, y_n = y0, px_n = px0, py_n = py0;
    double x_n_next, y_n_next, px_n_next, py_n_next;
    double x_err, y_err, px_err, py_err, err;
    double k1x, k1y, k1px, k1py, k2x, k2y, k2px, k2py, k3x, k3y, k3px, k3py, k4x, k4y, k4px, k4py, k5x, k5y, k5px, k5py, k6x, k6y, k6px, k6py, k7x, k7y, k7px, k7py;
    double h = 0.01, p = 5, facmin = 0.7, facmax = 1.3, fac = 0.98;
    int step = 0;
    // coefficients:
    double  c1 = 0, c2 = (double)1 / 5, c3 = (double)3 / 10, c4 = (double)4 / 5, c5 = (double)8 / 9, c6 = 1, c7 = 1;
    double a21 = (double)1 / 5;
    double a31 = (double)3 / 40, a32 = (double)9 / 40;
    double a41 = (double)44 / 45, a42 = (double)-56 / 15, a43 = (double)32 / 9;
    double a51 = (double)19372 / 6561, a52 = (double)-25360 / 2187, a53 = (double)64448 / 6561, a54 = (double)-212 / 729;
    double a61 = (double)9017 / 3168, a62 = (double)-355 / 33, a63 = (double)46732 / 5247, a64 = (double)49 / 176, a65 = (double)-5103 / 18656;
    double a71 = (double)35 / 384, a72 = (double)0, a73 = (double)500 / 1113, a74 = (double)125 / 192, a75 = (double)-2187 / 6784, a76 = (double)11 / 84;
    double  b1 = (double)35 / 384, b2 = (double)0, b3 = (double)500 / 1113, b4 = (double)125 / 192, b5 = (double)-2187 / 6784, b6 = (double)11 / 84, b7 = 0;
    double  d1 = (double)5179 / 57600, d2 = (double)0, d3 = (double)7571 / 16695, d4 = (double)393 / 640, d5 = (double)-92097 / 339200, d6 = (double)187 / 2100, d7 = (double)1 / 40;

    while (t_n < t_end) {
        step++;
        k1x = dot_x(t_n, x_n, y_n, px_n, py_n);
        k1y = dot_y(t_n, x_n, y_n, px_n, py_n);
        k1px = dot_px(t_n, x_n, y_n, px_n, py_n);
        k1py = dot_py(t_n, x_n, y_n, px_n, py_n);

        k2x = dot_x(t_n + c2 * h, x_n + h * (a21 * k1x), y_n + h * (a21 * k1y), px_n + h * (a21 * k1px), py_n + h * (a21 * k1py));
        k2y = dot_y(t_n + c2 * h, x_n + h * (a21 * k1x), y_n + h * (a21 * k1y), px_n + h * (a21 * k1px), py_n + h * (a21 * k1py));
        k2px = dot_px(t_n + c2 * h, x_n + h * (a21 * k1x), y_n + h * (a21 * k1y), px_n + h * (a21 * k1px), py_n + h * (a21 * k1py));
        k2py = dot_py(t_n + c2 * h, x_n + h * (a21 * k1x), y_n + h * (a21 * k1y), px_n + h * (a21 * k1px), py_n + h * (a21 * k1py));

        k3x = dot_x(t_n + c3 * h, x_n + h * (a31 * k1x + a32 * k2x), y_n + h * (a31 * k1y + a32 * k2y), px_n + h * (a31 * k1px + a32 * k2px), py_n + h * (a31 * k1py + a32 * k2py));
        k3y = dot_y(t_n + c3 * h, x_n + h * (a31 * k1x + a32 * k2x), y_n + h * (a31 * k1y + a32 * k2y), px_n + h * (a31 * k1px + a32 * k2px), py_n + h * (a31 * k1py + a32 * k2py));
        k3px = dot_px(t_n + c3 * h, x_n + h * (a31 * k1x + a32 * k2x), y_n + h * (a31 * k1y + a32 * k2y), px_n + h * (a31 * k1px + a32 * k2px), py_n + h * (a31 * k1py + a32 * k2py));
        k3py = dot_py(t_n + c3 * h, x_n + h * (a31 * k1x + a32 * k2x), y_n + h * (a31 * k1y + a32 * k2y), px_n + h * (a31 * k1px + a32 * k2px), py_n + h * (a31 * k1py + a32 * k2py));

        k4x = dot_x(t_n + c4 * h, x_n + h * (a41 * k1x + a42 * k2x + a43 * k3x), y_n + h * (a41 * k1y + a42 * k2y + a43 * k3y), px_n + h * (a41 * k1px + a42 * k2px + a43 * k3px), py_n + h * (a41 * k1py + a42 * k2py + a43 * k3py));
        k4y = dot_y(t_n + c4 * h, x_n + h * (a41 * k1x + a42 * k2x + a43 * k3x), y_n + h * (a41 * k1y + a42 * k2y + a43 * k3y), px_n + h * (a41 * k1px + a42 * k2px + a43 * k3px), py_n + h * (a41 * k1py + a42 * k2py + a43 * k3py));
        k4px = dot_px(t_n + c4 * h, x_n + h * (a41 * k1x + a42 * k2x + a43 * k3x), y_n + h * (a41 * k1y + a42 * k2y + a43 * k3y), px_n + h * (a41 * k1px + a42 * k2px + a43 * k3px), py_n + h * (a41 * k1py + a42 * k2py + a43 * k3py));
        k4py = dot_py(t_n + c4 * h, x_n + h * (a41 * k1x + a42 * k2x + a43 * k3x), y_n + h * (a41 * k1y + a42 * k2y + a43 * k3y), px_n + h * (a41 * k1px + a42 * k2px + a43 * k3px), py_n + h * (a41 * k1py + a42 * k2py + a43 * k3py));

        k5x = dot_x(t_n + c5 * h, x_n + h * (a51 * k1x + a52 * k2x + a53 * k3x + a54 * k4x), y_n + h * (a51 * k1y + a52 * k2y + a53 * k3y + a54 * k4y), px_n + h * (a51 * k1px + a52 * k2px + a53 * k3px + a54 * k4px), py_n + h * (a51 * k1py + a52 * k2py + a53 * k3py + a54 * k4py));
        k5y = dot_y(t_n + c5 * h, x_n + h * (a51 * k1x + a52 * k2x + a53 * k3x + a54 * k4x), y_n + h * (a51 * k1y + a52 * k2y + a53 * k3y + a54 * k4y), px_n + h * (a51 * k1px + a52 * k2px + a53 * k3px + a54 * k4px), py_n + h * (a51 * k1py + a52 * k2py + a53 * k3py + a54 * k4py));
        k5px = dot_px(t_n + c5 * h, x_n + h * (a51 * k1x + a52 * k2x + a53 * k3x + a54 * k4x), y_n + h * (a51 * k1y + a52 * k2y + a53 * k3y + a54 * k4y), px_n + h * (a51 * k1px + a52 * k2px + a53 * k3px + a54 * k4px), py_n + h * (a51 * k1py + a52 * k2py + a53 * k3py + a54 * k4py));
        k5py = dot_py(t_n + c5 * h, x_n + h * (a51 * k1x + a52 * k2x + a53 * k3x + a54 * k4x), y_n + h * (a51 * k1y + a52 * k2y + a53 * k3y + a54 * k4y), px_n + h * (a51 * k1px + a52 * k2px + a53 * k3px + a54 * k4px), py_n + h * (a51 * k1py + a52 * k2py + a53 * k3py + a54 * k4py));

        k6x = dot_x(t_n + c6 * h, x_n + h * (a61 * k1x + a62 * k2x + a63 * k3x + a64 * k4x + a65 * k5x), y_n + h * (a61 * k1y + a62 * k2y + a63 * k3y + a64 * k4y + a65 * k5y), px_n + h * (a61 * k1px + a62 * k2px + a63 * k3px + a64 * k4px + a65 * k5px), py_n + h * (a61 * k1py + a62 * k2py + a63 * k3py + a64 * k4py + a65 * k5py));
        k6y = dot_y(t_n + c6 * h, x_n + h * (a61 * k1x + a62 * k2x + a63 * k3x + a64 * k4x + a65 * k5x), y_n + h * (a61 * k1y + a62 * k2y + a63 * k3y + a64 * k4y + a65 * k5y), px_n + h * (a61 * k1px + a62 * k2px + a63 * k3px + a64 * k4px + a65 * k5px), py_n + h * (a61 * k1py + a62 * k2py + a63 * k3py + a64 * k4py + a65 * k5py));
        k6px = dot_px(t_n + c6 * h, x_n + h * (a61 * k1x + a62 * k2x + a63 * k3x + a64 * k4x + a65 * k5x), y_n + h * (a61 * k1y + a62 * k2y + a63 * k3y + a64 * k4y + a65 * k5y), px_n + h * (a61 * k1px + a62 * k2px + a63 * k3px + a64 * k4px + a65 * k5px), py_n + h * (a61 * k1py + a62 * k2py + a63 * k3py + a64 * k4py + a65 * k5py));
        k6py = dot_py(t_n + c6 * h, x_n + h * (a61 * k1x + a62 * k2x + a63 * k3x + a64 * k4x + a65 * k5x), y_n + h * (a61 * k1y + a62 * k2y + a63 * k3y + a64 * k4y + a65 * k5y), px_n + h * (a61 * k1px + a62 * k2px + a63 * k3px + a64 * k4px + a65 * k5px), py_n + h * (a61 * k1py + a62 * k2py + a63 * k3py + a64 * k4py + a65 * k5py));

        k7x = dot_x(t_n + c7 * h, x_n + h * (a71 * k1x + a72 * k2x + a73 * k3x + a74 * k4x + a75 * k5x + a76 * k6x), y_n + h * (a71 * k1y + a72 * k2y + a73 * k3y + a74 * k4y + a75 * k5y + a76 * k6y), px_n + h * (a71 * k1px + a72 * k2px + a73 * k3px + a74 * k4px + a75 * k5px + a76 * k6px), py_n + h * (a71 * k1py + a72 * k2py + a73 * k3py + a74 * k4py + a75 * k5py + a76 * k6py));
        k7y = dot_y(t_n + c7 * h, x_n + h * (a71 * k1x + a72 * k2x + a73 * k3x + a74 * k4x + a75 * k5x + a76 * k6x), y_n + h * (a71 * k1y + a72 * k2y + a73 * k3y + a74 * k4y + a75 * k5y + a76 * k6y), px_n + h * (a71 * k1px + a72 * k2px + a73 * k3px + a74 * k4px + a75 * k5px + a76 * k6px), py_n + h * (a71 * k1py + a72 * k2py + a73 * k3py + a74 * k4py + a75 * k5py + a76 * k6py));
        k7px = dot_px(t_n + c7 * h, x_n + h * (a71 * k1x + a72 * k2x + a73 * k3x + a74 * k4x + a75 * k5x + a76 * k6x), y_n + h * (a71 * k1y + a72 * k2y + a73 * k3y + a74 * k4y + a75 * k5y + a76 * k6y), px_n + h * (a71 * k1px + a72 * k2px + a73 * k3px + a74 * k4px + a75 * k5px + a76 * k6px), py_n + h * (a71 * k1py + a72 * k2py + a73 * k3py + a74 * k4py + a75 * k5py + a76 * k6py));
        k7py = dot_py(t_n + c7 * h, x_n + h * (a71 * k1x + a72 * k2x + a73 * k3x + a74 * k4x + a75 * k5x + a76 * k6x), y_n + h * (a71 * k1y + a72 * k2y + a73 * k3y + a74 * k4y + a75 * k5y + a76 * k6y), px_n + h * (a71 * k1px + a72 * k2px + a73 * k3px + a74 * k4px + a75 * k5px + a76 * k6px), py_n + h * (a71 * k1py + a72 * k2py + a73 * k3py + a74 * k4py + a75 * k5py + a76 * k6py));

        x_n_next = x_n + h * (b1 * k1x + b2 * k2x + b3 * k3x + b4 * k4x + b5 * k5x + b6 * k6x + b7 * k7x);
        x_err = x_n + h * (d1 * k1x + d2 * k2x + d3 * k3x + d4 * k4x + d5 * k5x + d6 * k6x + d7 * k7x);
        x_err = fabs(x_n_next - x_err);

        y_n_next = y_n + h * (b1 * k1y + b2 * k2y + b3 * k3y + b4 * k4y + b5 * k5y + b6 * k6y + b7 * k7y);
        y_err = y_n + h * (d1 * k1y + d2 * k2y + d3 * k3y + d4 * k4y + d5 * k5y + d6 * k6y + d7 * k7y);
        y_err = fabs(y_n_next - y_err);

        px_n_next = px_n + h * (b1 * k1px + b2 * k2px + b3 * k3px + b4 * k4px + b5 * k5px + b6 * k6px + b7 * k7px);
        px_err = px_n + h * (d1 * k1px + d2 * k2px + d3 * k3px + d4 * k4px + d5 * k5px + d6 * k6px + d7 * k7px);
        px_err = fabs(px_n_next - px_err);

        py_n_next = py_n + h * (b1 * k1py + b2 * k2py + b3 * k3py + b4 * k4py + b5 * k5py + b6 * k6py + b7 * k7py);
        py_err = py_n + h * (d1 * k1py + d2 * k2py + d3 * k3py + d4 * k4py + d5 * k5py + d6 * k6py + d7 * k7py);
        py_err = fabs(py_n_next - py_err);

        err = max(max(max(x_err, y_err), px_err), py_err);

        if (err <= tol) {
            x_n = x_n_next;
            y_n = y_n_next;
            px_n = px_n_next;
            py_n = py_n_next;
            t_n += h;
            h = h * min(facmax, max(facmin, fac * pow((tol / err), (double)(1 / (p + 1)))));
        }
        else {
            h = h * min(facmax, max(facmin, fac * pow((tol / err), (double)(1 / (p + 1)))));
        }
        if (h < tol || step > 10e5) {
            printf("\n\tToo many steps in shoot() function or h = %e is too small\n\n", h);
            exit(0);
        }
        if (t_n + h > t_end) {
            h = t_end - t_n;
        }
        //printf("x_n = %.2e\ty_n = %.2e\tpx_n = %.2e\tpy_n = %.2e\n", x_n, y_n, px_n, py_n);
    }
    ret_arr[0] = x_n;
    ret_arr[1] = y_n;
    ret_arr[2] = px_n;
    ret_arr[3] = py_n;
    //printf("\n\n\t\t\t x_n = %f\ty_n = %f\ty0 = %f\tpx0 = %f", ret_arr[0], ret_arr[1], y0, px0);
};
double dX1_da1(double x0, double y0, double alpha1, double alpha2, double t0, double T) {
    double initial_values[4], incremented_values[4];
    shoot(x0, y0, alpha1, alpha2, t0, T, initial_values);
    shoot(x0, y0, alpha1 + tol, alpha2, t0, T, incremented_values);
    double derivative = (initial_values[1] - incremented_values[1]) / tol; // (y(T) - y~(T)) / delta
    return derivative;
}
double dX1_da2(double x0, double y0, double alpha1, double alpha2, double t0, double T) {
    double initial_values[4], incremented_values[4];
    shoot(x0, y0, alpha1, alpha2, t0, T, initial_values);
    shoot(x0, y0, alpha1, alpha2 + tol, t0, T, incremented_values);
    double derivative = (initial_values[1] - incremented_values[1]) / tol; // (y(T) - y~(T)) / delta
    return derivative;
}
double dX2_da1(double x0, double y0, double alpha1, double alpha2, double t0, double T) {
    double initial_values[4], incremented_values[4];
    shoot(x0, y0, alpha1, alpha2, t0, T, initial_values);
    shoot(x0, y0, alpha1 + tol, alpha2, t0, T, incremented_values);
    double derivative = (initial_values[2] - incremented_values[2]) / tol; // (px(T) - px~(T)) / delta
    return derivative;
}
double dX2_da2(double x0, double y0, double alpha1, double alpha2, double t0, double T) {
    double initial_values[4], incremented_values[4];
    shoot(x0, y0, alpha1, alpha2, t0, T, initial_values);
    shoot(x0, y0, alpha1, alpha2 + tol, t0, T, incremented_values);
    double derivative = (initial_values[2] - incremented_values[2]) / tol; // (px(T) - px~(T)) / delta
    return derivative;
}

void print_solution_vector(double X[N]) {
    cout << "The solution for the system:" << endl;
    for (int i = 0; i < N; i++) {
        cout.precision(2);
        cout << fixed << "x" << i + 1 << " = " << X[i] << "\t";
        if (i % 15 == 0) {
            cout << endl;
        }
    }
    cout << endl;
}
void print_matrix(double A[N][N + 1]) {
    cout.precision(2);
    for (int i = 0; i < N; i++) {
        cout << fixed << endl;
        for (int j = 0; j <= N; j++) {
            cout << A[i][j] << "\t";
        }
    }
    cout << endl;
}
void swap_row(double A[N][N + 1], int i, int j) {
    for (int k = 0; k <= N; k++) {
        double temp = A[i][k];
        A[i][k] = A[j][k];
        A[j][k] = temp;
    }
}
int forward_elimination(double A[N][N + 1]) {
    for (int k = 0; k < N; k++) {
        int index_of_max_element_in_row = k;
        int value_of_max_element_in_row = fabs(A[index_of_max_element_in_row][k]);
        for (int i = k + 1; i < N; i++)
            if (fabs(A[i][k]) > value_of_max_element_in_row)
                value_of_max_element_in_row = A[i][k], index_of_max_element_in_row = i;
        if (fabs(A[k][index_of_max_element_in_row]) < EPSILON) {
            return k;
        }
        if (index_of_max_element_in_row != k) {
            swap_row(A, k, index_of_max_element_in_row);
        }
        for (int i = k + 1; i < N; i++) {
            double c = A[i][k] / A[k][k];
            for (int j = k + 1; j <= N; j++)
                A[i][j] -= c * A[k][j];
            A[i][k] = 0;
        }
    }
    return -1;
}
void back_substitution(double A[N][N + 1], double X[N]) {
    for (int i = N - 1; i >= 0; i--) {
        X[i] = A[i][N];
        for (int j = i + 1; j < N; j++) {
            X[i] -= A[i][j] * X[j];
        }
        X[i] = X[i] / A[i][i];
    }
}
void Gauss(double A[N][N + 1], double X[N]) {
    int ld_flag = forward_elimination(A);
    if (ld_flag != -1) {
        cout << "At least two rows of the matrix are linearly dependent";
        return;
    }
    back_substitution(A, X);
}

int main() {
    // given values from initial system
    double x0 = 0, y0 = 1;
    double t0 = 0, tn = 1;
    double yT = 0, pxT = 0;
    // initial parameters
    double alpha1 = 1, alpha2 = 1; // px(0), py(0)
    double initial_values[4], next_values[4]; // for a given iteration, NOT general
    // auxiliary variables:
    double dX1_dalpha1, dX1_dalpha2, dX2_dalpha1, dX2_dalpha2, a, b, c, d, jacobian[2][3], inverse_determinant, inverse_jacobian[2][2], hn[2];
    double alpha_next[2], f_next[2], X[2];
    int steps = 0, k;
    shoot(x0, y0, alpha1, alpha2, t0, tn, initial_values);
    double f_initial[2] = { yT - initial_values[1], pxT - initial_values[2] };
    while (sqrt(f_initial[0] * f_initial[0] + f_initial[1] * f_initial[1]) > tol) {
        steps++;
        if (steps > 10e3) {
            printf("\n\n\t\titeration limit exceeded allotted value, exiting program...\n\n");
            exit(0);
        }
        shoot(x0, y0, alpha1, alpha2, t0, tn, initial_values);
        //printf("\ninitial_values:\n\tx = %.2f\ty = %.2f", initial_values[0], initial_values[1]);
        f_initial[0] = yT - initial_values[1];
        f_initial[1] = pxT - initial_values[2];
        a = dX1_da1(x0, y0, alpha1, alpha2, t0, tn);
        b = dX1_da2(x0, y0, alpha1, alpha2, t0, tn);
        c = dX2_da1(x0, y0, alpha1, alpha2, t0, tn);
        d = dX2_da2(x0, y0, alpha1, alpha2, t0, tn);
        jacobian[0][0] = a;
        jacobian[0][1] = b;
        jacobian[1][0] = c;
        jacobian[1][1] = d;
        /*
        inverse_determinant = 1 / (a * d - b * c);
        inverse_jacobian[0][0] = d * inverse_determinant;
        inverse_jacobian[0][1] = -b * inverse_determinant;
        inverse_jacobian[1][0] = -c * inverse_determinant;
        inverse_jacobian[1][1] = a * inverse_determinant;
        hn[0] = inverse_jacobian[0][0] * f_initial[0] + inverse_jacobian[0][1] * f_initial[1];
        hn[1] = inverse_jacobian[1][0] * f_initial[0] + inverse_jacobian[1][1] * f_initial[1];
        */
        // jacobian * hn = f_initial - matrix-wise, find hn vector:
        // augmented matrix (jacobian|f_initial):
        jacobian[0][2] = f_initial[0];
        jacobian[1][2] = f_initial[1];
        //print_matrix(jacobian);
        Gauss(jacobian, hn);
        
        //print_solution_vector(hn);
        for (k = 0; k < 60; k++) {
            alpha_next[0] = alpha1 - pow(2., -k) * hn[0];
            alpha_next[1] = alpha2 - pow(2., -k) * hn[1];
            //printf("\nalpha_next1 = %f\talpha_next2 = %f\n", alpha_next[0], alpha_next[1]);
            shoot(x0, y0, alpha_next[0], alpha_next[1], t0, tn, next_values);
            //printf("\n***************************************************************************************************************\nnext_values:\tx = %.2f\ty = %.2f", next_values[0], next_values[1]);
            f_next[0] = yT - next_values[1];
            f_next[1] = pxT - next_values[2];
            //printf("\nsqrt_next = %e, sqrt_cur = %e", sqrt(f_next[0] * f_next[0] + f_next[1] * f_next[1]), sqrt(f_initial[0] * f_initial[0] + f_initial[1] * f_initial[1]));
            if (sqrt(f_next[0] * f_next[0] + f_next[1] * f_next[1]) < sqrt(f_initial[0] * f_initial[0] + f_initial[1] * f_initial[1])) {
                //if (fabs(f_next[0]) <= fabs(f_initial[0]) && fabs(f_next[1]) <= fabs(f_initial[1])) {
                alpha1 = alpha_next[0];
                alpha2 = alpha_next[1];
                //printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\t\tBREAKING! on k = %i\n\n\n", k);
                break;
            }
        }
        //printf("\nxT=%.2f\tyT=%.2f\t\twhile_sqrt=%e\twhile>tol:%i\n", f_next[0], f_next[1], sqrt(f_initial[0] * f_initial[0] + f_initial[1] * f_initial[1]), sqrt(f_initial[0] * f_initial[0] + f_initial[1] * f_initial[1]) > tol);

    };
    printf("\n\tFound solution for Boundry Value Problem (Alpha = %.1f) in %i iterations:\n\tWhen px(0) = %.17e\tand\tpy(0) = %.17e,\n\t     y(T) = %.17e\tand\tpx(T) = %.17e\n\n", beta, steps, alpha1, alpha2, next_values[1], next_values[2]);
    /*
    double A[N][N + 1] = {  {  3,   4,      -9,     5,      -14},   // N = 4, -8, -2, -2, 0
                            {-15,   -12,    50,     -16,    44},
                            {-27,   -36,    73,      8,     142},
                            {  9,   12,     -10,    -16,    -76} };
    double X[N];
    Gauss(A, X);
    for (int i = 0; i < N; i++) {
        printf("\n\tx[%i] = %f", i, X[i]);
    }
    */
}