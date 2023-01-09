#include "fnder2.h"
//#include "spmak2.h"

/*
function nsp = fnder2(sp, n)

*/
void  fnder2(stusp1* nsp, const stusp1* sp, int  n)
{
    /*
    knot = sp.knots;
    coef = sp.coefs;
    p = length(knot) - sp.number - 1;

    knotnumber = length(knot);
    coefnumber = length(coef(1, :));
    */

    Matrix  *knot, * coef;
    int i, p, knotnumber, coefnumber;

    i = 0;
    p = 0;
    knot = sp->knots;
    coef = sp->coefs;
    p          = knot->column - sp->number - 1;
    knotnumber = knot->column;
    coefnumber = coef->column;

    /*
    if n == 0
        nknot = knot;
    ncoef = coef;
    elseif n < p + 1
        nknot = knot(n + 1:knotnumber - n);
    for k = 1:n
        for i = 1 : coefnumber - k
            if knot(i + p + 1) == knot(i + k) % 处理重节点分母为0的情况
                c = 0;
            else
                c = (p + 1 - k) / (knot(i + p + 1) - knot(i + k));
    end
        P = coef(:, i + 1) - coef(:, i);
    ncoef(:, i) = c * P;
    end
        coef = ncoef(:, 1 : coefnumber - k);
    end
        ncoef = ncoef(:, 1 : coefnumber - n);
    else% 求导超过p + 1次，
        ncoef = zeros(length(coef(:, 1)), coefnumber - n); % 不报错，输出0样条
    end
    */

    Matrix  *nknot, *ncoef, *P;
    int k, j, m;
    double c;

    k = 0;
    c = 0;
    j = 0;
    m = 0;

    if ( n == 0 )
    {
        nknot = knot;
        ncoef = coef;
    }
    else if ( n < p + 1 )
    {
        P         = (Matrix*)malloc(sizeof(Matrix));
        P->data   = (double*)malloc(coef->row * sizeof(double));
        P->row    = 1;
        P->column = coef->row;
        nknot       = (Matrix*)malloc(sizeof(Matrix));
        nknot->data = (double*)malloc((knotnumber - n - (n + 1) + 1) * sizeof(double));
        //nknot = knot(n + 1:knotnumber - n);
        nknot->row    = knot->row;
        nknot->column = knotnumber - n - (n + 1)+1;
        for (i=n+1; i <= knotnumber-n; i++)
        {
            nknot->data[i-n-1] = knot->data[i-1];
        }
        for ( k = 1; k<= n; k++  )
        {
            ncoef         = (Matrix*)malloc(sizeof(Matrix));
            ncoef->data   = (double*)malloc(coef->row * (coefnumber - k) * sizeof(double));
            ncoef->row    = coef->row;
            ncoef->column = coefnumber - k;
            for ( i = 1; i <= coefnumber - k; i++ )
            {
                //if knot(i + p + 1) == knot(i + k)
                if ( (knot->data[(i + k - 1)] - 0.00005) <= knot->data[(i + p + 1 - 1)] <= knot->data[(i + k - 1)] + 0.00005 )
                {
                    c = 0;
                }
                else
                {
                    c = (p + 1 - k) / (knot->data[(i + p + 1 - 1)] - knot->data[(i + k - 1)] ) ;
                }
                //P = coef(:, i + 1) - coef(:, i);
                for ( j = 0; j < coef->row; j++ )
                {
                    P->data[j] = coef->data[((i+1-1) * coef->row)+ j] - coef->data[((i - 1) * coef->row) + j];
                }
                //ncoef(:, i) = c * P;
                for (j = 0; j < coef->row; j++)
                {
                    ncoef->data[(i-1) * coef->row + j] = c * P->data[j];
                }
            }

            //coef = ncoef(:, 1 : coefnumber - k);
            coef->row    = ncoef->row;
            coef->column = coefnumber - k;
            for (m = 0; m < coefnumber - k; m++)
            {
                for (j = 0; j < ncoef->row; j++)
                {
                    coef->data[ m * ncoef->row + j] = ncoef->data[ m * ncoef->row + j];
                }
            }

        }

        //ncoef = ncoef(:, 1 : coefnumber - n);
        ncoef->column = coefnumber - n;
        for (m = 0; m < coefnumber - n; m++)
        {
            for (j = 0; j < ncoef->row; j++)
            {
                ncoef->data[m * ncoef->row + j] = ncoef->data[m * ncoef->row + j];
            }

        }

    }
    else
    {
        //ncoef = zeros(length(coef(:, 1)), coefnumber - n); % 不报错，输出0样条
        nknot = (Matrix*)malloc(sizeof(Matrix));
        nknot->data    = (double*)malloc((knotnumber - n - (n + 1) + 1) * sizeof(double));
        nknot->row     = 1;
        nknot->column  = knotnumber - n - (n + 1) + 1;
        ncoef          = Matrix_gen1(coef->row, coefnumber - n, 0.00 );
    }

    //nsp = spmak2(nknot, ncoef);
    nsp->number = ncoef->column;
    nsp->knots  = nknot;
    nsp->coefs  = ncoef;
    //spmak2(nknot, ncoef, nsp);

    M_free(P);
}





