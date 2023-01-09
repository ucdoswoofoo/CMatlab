#include "Binsert.h"
#include "ZXMat.h"
#include "matrix.h"


void  Binsert(Matrix* ncontr, const stusp1* sp, const Matrix* superKnots)
{
    Matrix *ocontr, *oknots, * nknots;
    int k, n, i, j, m, secin, w, x, ii, jj;
    double alpha;

    k = 0;
    n = 0;
    i = 0;
    j = 0;
    m = 0;
    k = 0;
    w = 0;
    x = 0;
    ii = 0;
    jj = 0;
    secin = 0;
    alpha = 0.00;

    //ocontr = sp.coefs;% 旧的控制点
    //oknots = sp.knots;% 旧的节点
    ocontr = Matrix_copy(sp->coefs);
    oknots = Matrix_copy(sp->knots);

    //k = 3; % B样条曲线的次数
    k = 3;
    //n = sp.number;
    n = sp->number;
    //m = length(superKnots);
    m = superKnots->column;

    //ncontr = zeros(5, n + m);  需要再函数外部实现
    //ncontr = Matrix_gen1(5, n + m, 0.00);
    //nknots = Matrix_gen1(1, oknots->column + 1, 0.00);
    //secin = 0;
    //for j = 1:m
    nknots = NULL;
    for (j=1; j<=m; j++)
    {
        //for i = 4 : n + 1 + secin
        for (i=4; i<= n+1+secin; i++ )
        {
            //if superKnots(j) > oknots(i) && superKnots(j) < oknots(i+1)
            if ((superKnots->data[j-1] > oknots->data[i-1]) && (superKnots->data[j-1] < oknots->data[i + 1-1]))
            {
                //secin = secin + 1;
                secin = secin + 1;
                if ( nknots )
                {
                    jj = nknots->column;
                    M_free( nknots );
                }
                //nknots = Matrix_gen1(1, oknots->column + 1, 0.00);
                nknots = Matrix_gen1(1, i + 1 + oknots->column - i, 0.00);
                //nknots = [oknots(1:i) superKnots(j) oknots(i + 1:end)];
                for (ii=0; ii<i; ii++)
                {
                    nknots->data[ii] = oknots->data[ii];
                }
                //nknots->data[0] = oknots->data[(i-1) * oknots->row];
                nknots->data[i] = superKnots->data[j - 1];

                //nknots->data[2] = oknots->data[oknots->row * (oknots->column - 1) + i];
                for (ii = i; ii < oknots->column; ii++)
                {
                    nknots->data[ii+1] = oknots->data[ii];
                }

                //for w = 1:i - k
                for (w=1; w<=i-k; w++)
                {
                    //ncontr(:, w) = ocontr(:, w);
                    for (x = 0; x < ocontr->row; x++)
                    {
                        ncontr->data[(w-1)* ocontr->row + x] = ocontr->data[(w-1)* ocontr->row + x];
                    }
                }
                //for w = i - k + 1:i
                for (w = i-k+1; w <= i; w++)
                {
                    //alpha = (superKnots(j) - oknots(w)) / (oknots(w + k) - oknots(w));
                    alpha = (superKnots->data[j-1] - oknots->data[w-1]) / (oknots->data[w + k-1] - oknots->data[w-1]);
                    //ncontr(:, w) = (1 - alpha) * ocontr(:, w - 1) + alpha * ocontr(:, w);
                    for (x = 0; x < ocontr->row; x++)
                    {
                        ncontr->data[(w-1)* ncontr->row + x] = (1 - alpha) * ocontr->data[(w - 1 -1)* ocontr->row + x] + alpha * ocontr->data[(w-1)* ocontr->row + x];
                    }
                }
                //for w = i + 1:n + j
                for (w = i+1; w <= n + j; w++)
                {
                    //ncontr(:, w) = ocontr(:, w - 1);
                    for (x=0; x< ncontr->row; x++)
                    {
                        ncontr->data[(w - 1) * ncontr->row + x] = ocontr->data[(w-1-1) * ocontr->row + x];
                    }
                }
                M_free(ocontr);
                ocontr = Matrix_copy(ncontr);
                //oknots = nknots;
                M_free(oknots);
                oknots = Matrix_copy(nknots);
                M_free(nknots);
                nknots = NULL;
            }
        }
    }

   //% ncontr = ncontr(:, any(ncontr, 1));

}



