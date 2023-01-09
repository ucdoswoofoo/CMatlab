#include "fnval2.h"
#include "BaseFunction.h"

/*
% 输入：1 样条结构体f，由节点knots，控制点coefs，控制点个数number构成；2 参数数组s
% 输出：样条函数f在参数s处的取值

% 主函数
function v = fnval2(sp, s)
*/
void  fnval2(Matrix *v, const stusp1* sp, const double * s, int sLength)
{
    /*
    knot = sp.knots;
    coef = sp.coefs;
    n = sp.number;
    k = length(knot) - n - 1;

    m = length(coef(:, 1));
    r = length(s);
    v = zeros(m, r);
    */

    Matrix *knot, *coef;
    int n, m, k, i, ii, r, mm, nn;

    ii   = 0;
    knot = sp->knots;
    coef = sp->coefs;
    n    = sp->number;
    k    = knot->column - n - 1; 
    m    = coef->row;
    r    = sLength;

    for (i = 1; i <= r; i++)
    {
        //ii=find(knot<=s(i),1,'last'); %确定参数所在的节点区间
        for (mm = 0; mm < knot->column; mm++)
        {
            if (1 == knot->row)
            {
                //if ( ( (knot->data[mm] - 0.0005) <= s[i - 1]) && (s[i - 1] <= (knot->data[mm] + 0.0005) ))
                if ( knot->data[mm] <= s[i - 1] )
                {
                    ii = mm + 1;
                }
            }
            else
            {
                for (nn = 0; nn < knot->row; nn++)
                {
                    if (  knot->data[mm * knot->row + nn] <= s[i - 1] )
                    {
                        ii = mm * knot->row + nn + 1;
                    }
                }
            }
        }
        if ( 0 == ii )
        {
            ii = knot->row * knot->column;
        }
        if (ii == 1)
        {
            ii = k + 1;
        }


        if (ii == n + k + 1)
        {
            for (mm = 0; mm < coef->row; mm++)
            {
                v->data[ (i-1) * coef->row + mm] = coef->data[(coef->column - 1) * (coef->row) + mm];
            }
            continue;
        }

        if (k == 3)
        {
            //v(:, i) = coef(:, ii - 3) * BaseFunction(ii - 3, k, s(i), knot) + coef(:, ii - 2) * BaseFunction(ii - 2, k, s(i), knot) + coef(:, ii - 1) * BaseFunction(ii - 1, k, s(i), knot) + coef(:, ii) * BaseFunction(ii, k, s(i), knot);
            for (mm = 0; mm < coef->row; mm++)
            {
                v->data[coef->row * (i-1) + mm] = coef->data[(ii - 3-1) * coef->row + mm] * BaseFunction(ii - 3, k, s[i-1], knot->data ) + coef->data[coef->row * (ii - 2-1) + mm] * BaseFunction(ii - 2, k, s[(i-1)], knot->data) + coef->data[coef->row * (ii - 1-1) + mm] * BaseFunction(ii - 1, k, s[(i-1)], knot->data) + coef->data[coef->row * (ii-1) + mm] * BaseFunction(ii, k, s[(i-1)], knot->data);

            }

        }
        if (k == 2)
        {
            //v(:, i) = coef(:, ii - 2) * BaseFunction(ii - 2, k, s(i), knot) + coef(:, ii - 1) * BaseFunction(ii - 1, k, s(i), knot) + coef(:, ii) * BaseFunction(ii, k, s(i), knot);
            for (mm = 0; mm < coef->row; mm++)
            {
                v->data[coef->row * (i-1) + mm] = coef->data[coef->row * (ii - 2-1) + mm] * BaseFunction(ii - 2, k, s[(i-1)], knot->data) + coef->data[coef->row * (ii - 1-1) + mm] * BaseFunction(ii - 1, k, s[(i-1)], knot->data) + coef->data[coef->row * (ii-1) + mm] * BaseFunction(ii, k, s[(i-1)], knot->data);
            }
        }
        if (k == 1)
        {
            //v(:, i) = coef(:, ii - 1) * BaseFunction(ii - 1, k, s(i), knot) + coef(:, ii) * BaseFunction(ii, k, s(i), knot);
            for (mm = 0; mm < coef->row; mm++)
            {
                v->data[coef->row * (i-1) + mm] = coef->data[coef->row * (ii - 1-1) + mm] * BaseFunction(ii - 1, k, s[(i-1)], knot->data) + coef->data[coef->row * (ii-1) + mm] * BaseFunction(ii, k, s[i-1], knot->data);
            }
        }
        if (k == 0)
        {
            //v(:, i) = coef(:, ii) * BaseFunction(ii, k, s(i), knot);
            for (mm = 0; mm < coef->row; mm++)
            {
                v->data[coef->row * (i-1) + mm] = coef->data[coef->row * (ii-1) + mm] * BaseFunction(ii, k, s[i-1], knot->data);
            }
        }
    }

}





