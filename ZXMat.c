#include "ZXMat.h"




int ZXUnique(double  PrimVal[], int nCnt)
{
    int nNewLen = 0;
    int j = 0;

    for (int i = 0, j = 0; i < nCnt && j < nCnt; i++)
    {
        while (j < nCnt && PrimVal[i] == PrimVal[j])
        {
            j++;
        }

        if (j > i + 1 && j < nCnt)
        {
            PrimVal[i + 1] = PrimVal[j];
        }

        nNewLen++;
    }

    return nNewLen;
}



void ZXSort(double ArrVal[], int nCnt)
{
    int i, j;
    double t;

    for (i = nCnt - 2; i >= 0; i--)
    {
        for (j = 0; j < i; j++)
        {
            if (ArrVal[j] > ArrVal[j + 1])
            {
                t = ArrVal[j];
                ArrVal[j] = ArrVal[j + 1];
                ArrVal[j + 1] = t;
            }
        }
    }
}

void ZXdiff(double  PrimVal[], double DiffVal[], int nCnt)
{
    int i;

    i = 0;
    for (i = 0; i < nCnt - 1; i++)
    {
        DiffVal[i] = PrimVal[i + 1] - PrimVal[i];
    }
}

int ZXIsEmptyMat(Matrix _Mat)
{
    int IsEmpty, i;

    IsEmpty = 0;
    i = 0;

    if ( (_Mat.column = 0)  && (_Mat.row == 0) )
    {
        IsEmpty = 1;
    }

    return IsEmpty;

}



/*
// 得到矩阵的行列式
double ZXDeterminant(double** a, int k)
{
    double s = 1, det = 0, ** b;
    int i, j, m, n, c;


    b = (double**)calloc(k, sizeof(*b));
    for (i = 0; i < k; i++)
    {
        b[i] = (double*)calloc(k, sizeof(**b));
    }


    if (k == 1)
    {
        return (a[0][0]);
    }
    else
    {
        det = 0;
        for (c = 0; c < k; c++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < k; i++)
            {
                for (j = 0; j < k; j++)
                {
                    b[i][j] = 0;
                    if (i != 0 && j != c)
                    {
                        b[m][n] = a[i][j];
                        if (n < (k - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (a[0][c] * ZXDeterminant(b, k - 1));
            s = -1 * s;
        }
    }

    return (det);
}


// 得到矩阵特征码
void ZXCofactor(double** num, double** inverse, int f)
{
    double** b, ** fac;
    int p, q, m, n, i, j;
    //double** inverse;

    b = (double**)calloc(f, sizeof(*b));
    fac = (double**)calloc(f, sizeof(*fac));
    for (i = 0; i < f; i++)
    {
        b[i] = (double*)calloc(f, sizeof(**b));
        fac[i] = (double*)calloc(f, sizeof(**fac));
    }

    for (q = 0; q < f; q++)
    {
        for (p = 0; p < f; p++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < f; i++)
            {
                for (j = 0; j < f; j++)
                {
                    if (i != q && j != p)
                    {
                        b[m][n] = num[i][j];
                        if (n < (f - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1, q + p) * ZXDeterminant(b, f - 1);
        }
    }
    ZXTranspose(num, fac, inverse, f);
}


//得到逆矩阵
void ZXTranspose(double** num, double** fac, double** inverse, int r)
{
    int i, j;
    double** b,  d; //** inverse,

    b = (double**)calloc(r, sizeof(*b));
    //inverse = (double**)calloc(r, sizeof(*inverse));
    for (i = 0; i < r; i++)
    {
        b[i] = (double*)calloc(r, sizeof(**b));
        inverse[i] = (double*)calloc(r, sizeof(**inverse));
    }

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            b[i][j] = fac[j][i];
        }
    }

    d = ZXDeterminant(num, r);
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            inverse[i][j] = b[i][j] / d;
        }
    }

}

double ZXArrABSMax(double* Arr, int Count)
{
    return 0.0;
}


//矩阵与矩阵 相乘  要求Mat1的行数等于Mat2的列数
void ZXMatMUL(const int r1, const int c1, const CellMat m1, const int r2, const int c2, const CellMat m2, CellMat mul)
{
    int  i, j, k;

    if (r1 != c2)
    {
        return;
    }
    else 
    {
        for (i = 0; i < r1; i++)
        {
            for (j = 0; j < c2; j++)
            {
                mul.CellData[i][j] = 0;
                for (k = 0; k < c1; k++)
                {
                    mul.CellData[i][j] += m1.CellData[i][k] * m2.CellData[k][j];
                }
            }
        }
    }


}

//得到向量的欧几里得范数
double ZXMatnorm(CellMat Mat)
{
    int i;
    i = 0;
    double NormValue;

    NormValue = 0.00;

    for (i=0; i<Mat.iCol; i++)
    {
        NormValue = NormValue + pow(Mat.CellData[0][i], 2);
    }

    NormValue = pow(NormValue, 1 / 2);

    return NormValue;

}

CellMat ZXeye(int n)
{
    CellMat RetMat;
    int i, j;
    i = 0;
    j = 0;
    RetMat.CellData = (double**)calloc(n, sizeof(*RetMat.CellData));

    RetMat.iRow = n;
    RetMat.iCol = n;

    for (i = 0; i < n; i++)
    {
        RetMat.CellData[i] = (double*)calloc(RetMat.iCol, sizeof(**RetMat.CellData));
        for (j = 0; j < n; j++)
        {
            if ( i == j )
            {
                RetMat.CellData[i][j] = 1;
            }
            else
            {
                RetMat.CellData[i][j] = 0;
            }
        }
    }


    return RetMat;
}


double ZXDoubleMin(double* Arr, int Count)
{

    return 0.0;
}

void ZXABSMat(CellMat ZeroMat)
{
    int i, j;
    i = 0;
    j = 0;

    for (i = 0; i < ZeroMat.iRow; i++)
    {
        
        for (j = 0; j < ZeroMat.iCol; j++)
        {
            ZeroMat.CellData[i][j] = fabs(ZeroMat.CellData[i][j]);
        }
    }
}


//初始化矩阵:   新建矩阵(给出int Row, int Col)、赋初值InitValue
void  ZXMatInit(CellMat ZeroMat, int Row, int Col, double InitValue)
{
    int i, j;
    i = 0;
    j = 0;

    //ZeroMat.CellData = (double**)calloc(Row, sizeof(*ZeroMat.CellData));

    ZeroMat.iRow = Row;
    ZeroMat.iCol = Col;

    for (i=0; i< Row; i++)
    {
        ZeroMat.CellData[i] = (double*)calloc(ZeroMat.iCol, sizeof(**ZeroMat.CellData));
        for (j=0; j< Col; j++)
        {
            ZeroMat.CellData[i][j] = InitValue;
        }
    }
    //return ZeroMat.CellData;
}

//数字与矩阵相乘
void ZXMatSimMul(CellMat RetMat, CellMat ZeroMat, double MulValue)
{

}


//矩阵与矩阵 相加、减 Order: 0加 1减 
void ZXMatOper(CellMat RetMat, const CellMat Mat1, const CellMat Mat2, const int Order)
{

}

//截取矩阵的部分数据得到新的矩阵
void ZXGetMatPart(CellMat RetMat, const CellMat OrigMat, const int BeginRow, const int EndRow, const int BeginCol, const int EndCol)
{


}


//求矩阵的逆矩阵， 矩阵需要是方形矩阵
void ZXMatInv(CellMat RetMat, const CellMat PrimMat)
{
    double d;
    int n;
    n = PrimMat.iRow;
    d = ZXDeterminant(PrimMat.CellData, n);
    if (d == 0)
    {
        printf("Since the determinant is zerp (0), therefor inverse is not possible.");
    }
    else
    {
        ZXCofactor(PrimMat.CellData, RetMat.CellData, n);
    }
        
}
*/






