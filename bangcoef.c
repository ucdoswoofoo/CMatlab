#include "bangcoef.h"
#include "fnval2.h"

/*
% 根据bangbang控制求每个节点区间的放缩系数
% 输入：1 sp 样条结构体    2 Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax  1 * 1double型
% 输出：tknot新节点向量，talpha放缩系数
*/
void  bangcoef(Matrix* tknot, Matrix* talpha, Matrix* ttemp, const stusp1* sp, double Vmax, double Amax, double Jmax, double RVmax_A, double RVmax_C, double RAmax, double RJmax)
{
    /*
    n = sp.number;
    knot = sp.knots;
    dsp = fnder2(sp, 1);
    ddsp = fnder2(sp, 2);
    d3sp = fnder2(sp, 3);

    tknot = knot;
    talpha = 2 * ones(1, n - 3);
    ttemp = zeros(1, n - 3);
    */

    int n, i;
    Matrix *knot;
    stusp1* dsp, * ddsp, * d3sp;
    n    = 0;
    i    = 0;
    n    = sp->number;
    knot = sp->knots;

    dsp  = (stusp1*)malloc(sizeof(stusp1));
    ddsp = (stusp1*)malloc(sizeof(stusp1));
    d3sp = (stusp1*)malloc(sizeof(stusp1));

    //fnder2(stusp * nsp, const stusp * sp, int  n)
    stusp1 *sp1, *sp2, *sp3;

    sp1 = (stusp1*)malloc(sizeof(stusp1));
    sp2 = (stusp1*)malloc(sizeof(stusp1));
    sp3 = (stusp1*)malloc(sizeof(stusp1));

    sp1->number        = sp->number;
    sp1->coefs         = (Matrix*)malloc(sizeof(Matrix));
    sp1->coefs->column = sp->coefs->column;
    sp1->coefs->row    = sp->coefs->row;
    sp1->coefs->data   = (double*)malloc((sp1->coefs->column * sp1->coefs->row) * sizeof(double));
    for (i=0; i < (sp1->coefs->column * sp1->coefs->row); i++)
    {
        sp1->coefs->data[i] = sp->coefs->data[i];
    }
    sp1->knots         = (Matrix*)malloc(sizeof(Matrix));
    sp1->knots->column = sp->knots->column;
    sp1->knots->row    = sp->knots->row;
    sp1->knots->data   = (double*)malloc((sp1->knots->column * sp1->knots->row) * sizeof(double));
    for (i = 0; i < (sp1->knots->column * sp1->knots->row); i++)
    {
        sp1->knots->data[i] = sp->knots->data[i];
    }

    sp2->number        = sp->number;
    sp2->coefs         = (Matrix*)malloc(sizeof(Matrix));
    sp2->coefs->column = sp->coefs->column;
    sp2->coefs->row    = sp->coefs->row;
    sp2->coefs->data   = (double*)malloc((sp2->coefs->column * sp2->coefs->row) * sizeof(double));
    for (i = 0; i < (sp2->coefs->column * sp2->coefs->row); i++)
    {
        sp2->coefs->data[i] = sp->coefs->data[i];
    }
    sp2->knots         = (Matrix*)malloc(sizeof(Matrix));
    sp2->knots->column = sp->knots->column;
    sp2->knots->row    = sp->knots->row;
    sp2->knots->data   = (double*)malloc((sp2->knots->column * sp2->knots->row) * sizeof(double));
    for (i = 0; i < (sp2->knots->column * sp2->knots->row); i++)
    {
        sp2->knots->data[i] = sp->knots->data[i];
    }

    sp3->number        = sp->number;
    sp3->coefs         = (Matrix*)malloc(sizeof(Matrix));
    sp3->coefs->column = sp->coefs->column;
    sp3->coefs->row    = sp->coefs->row;
    sp3->coefs->data   = (double*)malloc((sp3->coefs->column * sp3->coefs->row) * sizeof(double));
    for (i = 0; i < (sp3->coefs->column * sp3->coefs->row); i++)
    {
        sp3->coefs->data[i] = sp->coefs->data[i];
    }
    sp3->knots         = (Matrix*)malloc(sizeof(Matrix));
    sp3->knots->column = sp->knots->column;
    sp3->knots->row    = sp->knots->row;
    sp3->knots->data   = (double*)malloc((sp3->knots->column * sp3->knots->row) * sizeof(double));
    for (i = 0; i < (sp3->knots->column * sp3->knots->row); i++)
    {
        sp3->knots->data[i] = sp->knots->data[i];
    }


    fnder2(dsp, sp1, 1);
    fnder2(ddsp, sp2, 2);
    fnder2(d3sp, sp3, 3);
    //fnder2(dsp, sp, 1);
    //fnder2(ddsp, sp, 2);
    //fnder2(d3sp, sp, 3);

    tknot = knot;

     for ( i = 4; i <= n; i++)
    {
        //t = knot(i) : (knot(i + 1) - knot(i)) / 100 : knot(i + 1);
        int j, m;
        Matrix t;

        j = 0;
        m = 1;
        for (j = 0; j < 35568; j++)
        {
            if ( knot->data[i-1] + ((knot->data[i] - knot->data[i-1]) / 100) * (m ) > knot->data[i] )
            {
                break;
            }
            m++;
        }
        

        t.data   = (double*)malloc((m)* sizeof(double));
        t.row    = 1;
        t.column = m;
        for (j = 0; j < m; j++)
        {
            t.data[j] = knot->data[i-1] + ((knot->data[i] - knot->data[i-1]) / 100) * (j);
        }

        //v = fnval2(dsp, t);
        Matrix *v;

        //v = (Matrix*)malloc(sizeof(Matrix));
        v = Matrix_gen1(dsp->coefs->row, t.column, 0.00);
        
        fnval2(v, dsp, t.data, t.column);
       
        double vx, vy, vz, va, vc;
        vx = 0.00;
        vy = 0.00;
        vz = 0.00;
        va = 0.00;
        vc = 0.00;
        /*
        vx = max(abs(v(1, :)));
        vy = max(abs(v(2, :)));
        vz = max(abs(v(3, :)));
        va = max(abs(v(4, :)));
        vc = max(abs(v(5, :)));
        */
        for (j = 0; j < v->column; j++)
        {
            if (vx < fabs(v->data[j * v->row ]))
            {
                vx = fabs(v->data[j * v->row]);
            }
        }
        for (j = 0; j < v->column; j++)
        {
            if (vy < fabs(v->data[ v->row * j + 1]))
            {
                vy = fabs(v->data[v->row * j + 1]);
            }
        }
        for (j = 0; j < v->column; j++)
        {
            if (vz < fabs(v->data[v->row * j + 2]))
            {
                vz = fabs(v->data[v->row * j + 2]);
            }
        }
        for (j = 0; j < v->column; j++)
        {
            if (va < fabs(v->data[v->row * j + 3]))
            {
                va = fabs(v->data[v->row * j + 3]);
            }
        }
        for (j = 0; j < v->column; j++)
        {
            if (vc < fabs(v->data[v->row * j  + 4]))
            {
                vc = fabs(v->data[v->row * j  + 4]);
            }
        }

        //alpha0 = min([Vmax / vx Vmax / vy Vmax / vz RVmax_A / va RVmax_C / vc]);
        Matrix Arr;
        double alpha0;

        alpha0 = 35568.00;
        Arr.row = 1;
        Arr.column = 5;
        Arr.data = (double*)malloc( Arr.column * sizeof( MATRIX_TYPE ) );

        Arr.data[0] = Vmax / vx;
        Arr.data[1] = Vmax / vy;
        Arr.data[2] = Vmax / vz;
        Arr.data[3] = RVmax_A / va;
        Arr.data[4] = RVmax_C / vc;

        for ( j = 0; j < 5; j++ )
        {
            if (alpha0 > Arr.data[j])
            {
                alpha0 = Arr.data[j];
            }
        }
        
        //double ataui, ataui1;
        Matrix  *atauiMat, *ataui1Mat, * atauiMatABS, * ataui1MatABS;
        //ataui = abs(fnval2(ddsp, knot(i)));
        //ataui1 = abs(fnval2(ddsp, knot(i + 1)));
        atauiMat  = Matrix_gen1(ddsp->coefs->row, 1, 0.00);
        ataui1Mat = Matrix_gen1(ddsp->coefs->row, 1, 0.00);

        fnval2(atauiMat, ddsp, &knot->data[i-1], 1);
        fnval2(ataui1Mat, ddsp, &knot->data[i], 1);
        
        atauiMatABS  = M_abs( atauiMat );
        ataui1MatABS = M_abs( ataui1Mat );
        
        //amin = min([sqrt(Amax / ataui(1)) sqrt(Amax / ataui(2)) sqrt(Amax / ataui(3)) sqrt(Amax / ataui1(1)) sqrt(Amax / ataui1(2)) sqrt(Amax / ataui1(3))...
        //    sqrt(RAmax / ataui(4)) sqrt(RAmax / ataui(5)) sqrt(RAmax / ataui1(4)) sqrt(RAmax / ataui1(5))]);

        Matrix Arr2;
        double amin;

        amin        = 35568.00;
        Arr2.row    = 1;
        Arr2.column = 10;
        Arr2.data   = (double*)calloc(Arr2.column, sizeof(MATRIX_TYPE));

        Arr2.data[0] = sqrt(Amax / atauiMatABS->data[0]);
        Arr2.data[1] = sqrt(Amax / atauiMatABS->data[1]);
        Arr2.data[2] = sqrt(Amax / atauiMatABS->data[2]);
        Arr2.data[3] = sqrt(Amax / ataui1MatABS->data[0]);
        Arr2.data[4] = sqrt(Amax / ataui1MatABS->data[1]);
        Arr2.data[5] = sqrt(Amax / ataui1MatABS->data[2]);
        Arr2.data[6] = sqrt(RAmax / atauiMatABS->data[3]);
        Arr2.data[7] = sqrt(RAmax / atauiMatABS->data[4]);
        Arr2.data[8] = sqrt(RAmax / ataui1MatABS->data[3]);
        Arr2.data[9] = sqrt(RAmax / ataui1MatABS->data[4]);

        for (j = 0; j < Arr2.column; j++)
        {
            if (amin > Arr2.data[j])
            {
                amin = Arr2.data[j];
            }
        }

        alpha0 = fmin(alpha0, amin);

        //jtaui = abs(fnval2(d3sp, knot(i + 1) / 2 + knot(i) / 2));
        Matrix  *jtauiMat, * jtauiMatABS, *fanval2s;
        fanval2s          = Matrix_gen1(1, 1, 0.00);
        fanval2s->data[0] = knot->data[i] / 2 + knot->data[i-1] / 2;
        jtauiMat = Matrix_gen1(d3sp->coefs->row, 1, 0.00);
        fnval2(jtauiMat, d3sp, &fanval2s->data[0], 1);
        jtauiMatABS = M_abs( jtauiMat );

        //jmin = min([(Jmax / jtaui(1)) ^ (1 / 3) (Jmax / jtaui(2)) ^ (1 / 3) (Jmax / jtaui(3)) ^ (1 / 3)...
        //    (RJmax / jtaui(4)) ^ (1 / 3) (RJmax / jtaui(5)) ^ (1 / 3)]);
        double jmin;
        Matrix Arr3;

        jmin = 35568.00;

        Arr3.row    = 1;
        Arr3.column = 5;
        Arr3.data   = (double*)calloc(Arr3.column, sizeof(MATRIX_TYPE));

        Arr3.data[0] = pow((Jmax / jtauiMatABS->data[0]), (1.0 / 3));
        Arr3.data[1] = pow((Jmax / jtauiMatABS->data[1]), (1.0 / 3));
        Arr3.data[2] = pow((Jmax / jtauiMatABS->data[2]), (1.0 / 3));
        Arr3.data[3] = pow((RJmax / jtauiMatABS->data[3]), (1.0 / 3));
        Arr3.data[4] = pow((RJmax / jtauiMatABS->data[4]), (1.0 / 3));

        //jmin = ZXDoubleMin(Arr3.CellData[0], 5);
        for (j = 0; j < Arr3.column; j++)
        {
            if (jmin > Arr3.data[j])
            {
                jmin = Arr3.data[j];
            }
        }
        
        alpha0 = fmin(alpha0, jmin);

        //talpha(i - 3) = alpha0;
        talpha->data[i - 3 - 1] = alpha0;

        if (i > 4)
        {
            //tknot(i) = (knot(i) - ttemp(i - 4)) / talpha(i - 4);
            //ttemp(i - 3) = knot(i) - talpha(i - 3) * tknot(i);
            tknot->data[i-1] = (knot->data[i-1] - ttemp->data[i - 4 - 1]) / talpha->data[i - 4 - 1];
            ttemp->data[i - 3 - 1] = knot->data[i - 1] - talpha->data[i - 3 - 1] * tknot->data[i - 1];
        }
       
        free(Arr3.data);
        free(Arr2.data);
        free(Arr.data);
        M_free(jtauiMat);
        M_free(jtauiMatABS);
        M_free(fanval2s);
        M_free(atauiMat);
        M_free(ataui1Mat);
        M_free(atauiMatABS);
        M_free(ataui1MatABS);
        M_free(v);

        free(sp1->coefs);
        free(sp1->knots);
        free(sp1);
        free(sp2->coefs);
        free(sp2->knots);
        free(sp2);
        free(sp3->coefs);
        free(sp3->knots);
        free(sp3);
        free(dsp);
        free(ddsp);
        free(d3sp);
        free(t.data);
        

    }
    //tknot(n + 1) = (knot(n + 1) - ttemp(n - 3)) / talpha(n - 3);
    //tknot(n + 2:end) = tknot(n + 1);
    tknot->data[(n + 1)-1] = (knot->data[(n + 1)-1] - ttemp->data[(n - 3)-1]) / talpha->data[(n - 3)-1];
    
    for (int i = n + 2; i <= tknot->column; i++)
    {
        tknot->data[i-1] = tknot->data[(n + 1)-1];
    }

    
}


