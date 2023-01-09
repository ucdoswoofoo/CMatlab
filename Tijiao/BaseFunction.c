#include "BaseFunction.h"


/* Function Definitions */
/*
% �����õĺ���
function Nik_u = BaseFunction(i, k, u, U)
% ���������Ni_k(u)��
% u������i�ڵ���ţ�k������������U�ڵ�ʸ��

* ���������Ni_k(u)��
 * u������i�ڵ���ţ�k������������U�ڵ�ʸ��
 *
 * Arguments    : double i
 *                double k
 *                double u
 *                const double U[8]
 * Return Type  : double
 */
double BaseFunction(double i, double k, double u, const double* U)
{
    double L1;
    double L2;
    double Nik_u;
    double alpha;
    double beta;
    /* �����õĺ��� */
    if (k == 0.0) {
        /* 0��B���� */
        if ((u >= U[(int)i - 1]) && (u < U[(int)i])) {
            Nik_u = 1.0;
        }
        else {
            Nik_u = 0.0;
        }
    }
    else {
        beta = i + k;
        alpha = U[(int)i - 1];
        L1 = U[(int)beta - 1] - alpha;
        /* ֧������ĳ��� */
        beta = U[(int)(beta + 1.0) - 1];
        L2 = beta - U[(int)i];
        if (L1 == 0.0) {
            /* �涨0/0=0 */
            alpha = 0.0;
        }
        else {
            alpha = (u - alpha) / L1;
        }
        if (L2 == 0.0) {
            beta = 0.0;
        }
        else {
            beta = (beta - u) / L2;
        }
        Nik_u = alpha * BaseFunction(i, k - 1.0, u, U) +
            beta * BaseFunction(i + 1.0, k - 1.0, u, U);
    }
    return Nik_u;
}

