
#include "globaltrans.h"

#include "bangc.h"
#include "bangcoef.h"
#include "fnval2.h"
#include "multal.h"

/*
% 整体调整曲线节点，使得曲线形状不变，时间延长或缩短, 曲线满足运动约束
% 输入：1 sp 样条结构体 2 Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax double型1 * 1
% 输出：nsp 样条结构体
function nsp = globaltrans(sp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax)
*/
void globaltrans(stusp1* nsp, const stusp1* sp, const double  Vmax, const double  Amax, const double  Jmax, const double RVmax_A, const double RVmax_C, const double  RAmax, const double  RJmax)
{
    //输入数据打印
    int iMat, yMat;
    iMat = 0;
    yMat = 0;
    printf("Vmax=%lf, Amax=%lf, Jmax=%lf, RVmax_A=%lf, RVmax_C=%lf, RAmax=%lf,  RJmax=%lf\n", Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);

    printf("sp->number = %d\n", sp->number);
    for (iMat=0; iMat < sp->coefs->column; iMat++)
    {
        for (yMat=0; yMat < sp->coefs->row ; yMat++)
        {
            printf("coefs Row=%d, Column=%d, Value=%lf\n", yMat, iMat, sp->coefs->data[sp->coefs->row * iMat + yMat]);
        }
    }

    for (iMat = 0; iMat < sp->knots->column ; iMat++)
    {
        for (yMat = 0; yMat < sp->knots->row ; yMat++)
        {
            printf("knots Row=%d, Column=%d, Value=%lf\n", yMat, iMat, sp->knots->data[sp->knots->row * iMat + yMat]);
        }
    }
    //输入数据打印
    //[~, alpha1, ~] = bangcoef(sp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);
    int i;
    Matrix * tknot1, *talpha, *ttemp, * tknot;

    i = 0;
    tknot1 = (Matrix*)malloc(sizeof(Matrix));
    talpha = (Matrix*)malloc(sizeof(Matrix));
    ttemp  = (Matrix*)malloc(sizeof(Matrix));

    //talpha = 2 * ones(1, n - 3);
    talpha  = Matrix_gen1(1, sp->number - 3, 1*2);
    //ttemp = zeros(1, n - 3);
    ttemp   = Matrix_gen1(1, sp->number - 3, 0);
    stusp1* sp1, * sp2, * sp3, * sp4, * sp5, * sp6; 

    sp1 = (stusp1*)malloc(sizeof(stusp1));
    sp2 = (stusp1*)malloc(sizeof(stusp1));
    sp3 = (stusp1*)malloc(sizeof(stusp1));
    sp4 = (stusp1*)malloc(sizeof(stusp1));
    sp5 = (stusp1*)malloc(sizeof(stusp1));
    sp6 = (stusp1*)malloc(sizeof(stusp1));

    sp1->number = sp->number;
    sp1->coefs = (Matrix*)malloc(sizeof(Matrix));
    sp1->coefs->column = sp->coefs->column;
    sp1->coefs->row = sp->coefs->row;
    sp1->coefs->data = (double*)malloc((sp1->coefs->column * sp1->coefs->row) * sizeof(double));
    for (i = 0; i < (sp1->coefs->column * sp1->coefs->row); i++)
    {
        sp1->coefs->data[i] = sp->coefs->data[i];
    }
    sp1->knots = (Matrix*)malloc(sizeof(Matrix));
    sp1->knots->column = sp->knots->column;
    sp1->knots->row = sp->knots->row;
    sp1->knots->data = (double*)malloc((sp1->knots->column * sp1->knots->row) * sizeof(double));
    for (i = 0; i < (sp1->knots->column * sp1->knots->row); i++)
    {
        sp1->knots->data[i] = sp->knots->data[i];
    }

    sp2->number = sp->number;
    sp2->coefs = (Matrix*)malloc(sizeof(Matrix));
    sp2->coefs->column = sp->coefs->column;
    sp2->coefs->row = sp->coefs->row;
    sp2->coefs->data = (double*)malloc((sp2->coefs->column * sp2->coefs->row) * sizeof(double));
    for (i = 0; i < (sp2->coefs->column * sp2->coefs->row); i++)
    {
        sp2->coefs->data[i] = sp->coefs->data[i];
    }
    sp2->knots = (Matrix*)malloc(sizeof(Matrix));
    sp2->knots->column = sp->knots->column;
    sp2->knots->row = sp->knots->row;
    sp2->knots->data = (double*)malloc((sp2->knots->column * sp2->knots->row) * sizeof(double));
    for (i = 0; i < (sp2->knots->column * sp2->knots->row); i++)
    {
        sp2->knots->data[i] = sp->knots->data[i];
    }

    sp3->number = sp->number;
    sp3->coefs = (Matrix*)malloc(sizeof(Matrix));
    sp3->coefs->column = sp->coefs->column;
    sp3->coefs->row = sp->coefs->row;
    sp3->coefs->data = (double*)malloc((sp3->coefs->column * sp3->coefs->row) * sizeof(double));
    for (i = 0; i < (sp3->coefs->column * sp3->coefs->row); i++)
    {
        sp3->coefs->data[i] = sp->coefs->data[i];
    }
    sp3->knots = (Matrix*)malloc(sizeof(Matrix));
    sp3->knots->column = sp->knots->column;
    sp3->knots->row = sp->knots->row;
    sp3->knots->data = (double*)malloc((sp3->knots->column * sp3->knots->row) * sizeof(double));
    for (i = 0; i < (sp3->knots->column * sp3->knots->row); i++)
    {
        sp3->knots->data[i] = sp->knots->data[i];
    }

    sp4->number = sp->number;
    sp4->coefs = (Matrix*)malloc(sizeof(Matrix));
    sp4->coefs->column = sp->coefs->column;
    sp4->coefs->row = sp->coefs->row;
    sp4->coefs->data = (double*)malloc((sp4->coefs->column * sp4->coefs->row) * sizeof(double));
    for (i = 0; i < (sp4->coefs->column * sp4->coefs->row); i++)
    {
        sp4->coefs->data[i] = sp->coefs->data[i];
    }
    sp4->knots = (Matrix*)malloc(sizeof(Matrix));
    sp4->knots->column = sp->knots->column;
    sp4->knots->row = sp->knots->row;
    sp4->knots->data = (double*)malloc((sp4->knots->column * sp4->knots->row) * sizeof(double));
    for (i = 0; i < (sp4->knots->column * sp4->knots->row); i++)
    {
        sp4->knots->data[i] = sp->knots->data[i];
    }

    sp5->number = sp->number;
    sp5->coefs = (Matrix*)malloc(sizeof(Matrix));
    sp5->coefs->column = sp->coefs->column;
    sp5->coefs->row = sp->coefs->row;
    sp5->coefs->data = (double*)malloc((sp5->coefs->column * sp5->coefs->row) * sizeof(double));
    for (i = 0; i < (sp5->coefs->column * sp5->coefs->row); i++)
    {
        sp5->coefs->data[i] = sp->coefs->data[i];
    }
    sp5->knots = (Matrix*)malloc(sizeof(Matrix));
    sp5->knots->column = sp->knots->column;
    sp5->knots->row = sp->knots->row;
    sp5->knots->data = (double*)malloc((sp5->knots->column * sp5->knots->row) * sizeof(double));
    for (i = 0; i < (sp5->knots->column * sp5->knots->row); i++)
    {
        sp5->knots->data[i] = sp->knots->data[i];
    }

    sp6->number = sp->number;
    sp6->coefs = (Matrix*)malloc(sizeof(Matrix));
    sp6->coefs->column = sp->coefs->column;
    sp6->coefs->row = sp->coefs->row;
    sp6->coefs->data = (double*)malloc((sp6->coefs->column * sp6->coefs->row) * sizeof(double));
    for (i = 0; i < (sp6->coefs->column * sp6->coefs->row); i++)
    {
        sp6->coefs->data[i] = sp->coefs->data[i];
    }
    sp6->knots = (Matrix*)malloc(sizeof(Matrix));
    sp6->knots->column = sp->knots->column;
    sp6->knots->row = sp->knots->row;
    sp6->knots->data = (double*)malloc((sp6->knots->column * sp6->knots->row) * sizeof(double));
    for (i = 0; i < (sp6->knots->column * sp6->knots->row); i++)
    {
        sp6->knots->data[i] = sp->knots->data[i];
    }

    bangcoef(tknot1, talpha, ttemp, sp1, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);

    //[tknot, alpha, beta] = multal(min(alpha1), sp);
    double alpha;
    Matrix * nalpha, *nbeta, * ncoef;
    //alpha = min(talpha);
    
    i = 0;
    alpha = talpha->data[0];
    for (i=0; i< talpha->column; i++)
    {
        if (alpha > talpha->data[i])
        {
            alpha = talpha->data[i];
        }
    }

    nalpha = Matrix_gen1(1, sp->number - 3, 1 * alpha);
    nbeta  = Matrix_gen1(1, sp->number - 3, 0);
    tknot  = Matrix_gen1(1, sp->number + 4, 0);

    multal(tknot,  nalpha, nbeta, alpha, sp2);
    //ncoef = bangc(sp, tknot, alpha, beta);
    ncoef = Matrix_gen1(5, sp->number, 0);
    bangc(ncoef, *sp3, *tknot, *nalpha, *nbeta);
    //nsp = spmak2(tknot, ncoef);
    //spmak2(tknot, ncoef, nsp);
    nsp->number = ncoef->column;
    nsp->knots = tknot;
    nsp->coefs = ncoef;

    //gap = sp.knots(end) / 10000;
    double gap;
    gap = sp->knots->data[sp->knots->column-1] / 10000;

    //t = 0:gap:sp.knots(end);
    //int i;
    i = 0;
    Matrix  t;

    //t = (Matrix)malloc(sizeof(Matrix));
    t.data   = (double*)malloc(10001 * sizeof(double));
    t.row    = 1;
    t.column = 10001;
    for (i = 0; i <= 10000; i++)
    {
        t.data[i] = gap  * i;
    }

    stusp1 * dsp, * ddsp, * d3sp;
    Matrix *v, *a, *jj;
 
    dsp  = (stusp1*)malloc(sizeof(stusp1));
    ddsp = (stusp1*)malloc(sizeof(stusp1));
    d3sp = (stusp1*)malloc(sizeof(stusp1));

    //% 速度
    // dsp = fnder2(sp, 1);
    fnder2( dsp, sp4, 1);
    //v = fnval2(dsp, t);
    v   = Matrix_gen1(sp->coefs->row, t.column, 0.00);
    fnval2(  v, dsp, t.data, t.column);

    //% 加速度
    // ddsp = fnder2(sp, 2);
    fnder2(ddsp, sp5, 2);
    //a = fnval2(ddsp, t);
    a   = Matrix_gen1(ddsp->coefs->row, t.column, 0.00);
    fnval2(a, ddsp, t.data, t.column);
    //% 加加速
    // d3sp = fnder2(sp, 3);
    fnder2(d3sp, sp6, 3);
    //jj = fnval2(d3sp, t);
    jj   = Matrix_gen1(d3sp->coefs->row, t.column, 0.00);
    fnval2(jj, d3sp, t.data, t.column);

    //alphav = max(abs([v(1, :) v(2, :) v(3, :)])) / Vmax;
    double alphav;
    alphav = 0.00;
    for (i = 0; i < v->column; i++)
    {
        if (alphav < fabs(v->data[v->row * i ]))
        {
            alphav = fabs(v->data[v->row * i ]);
        }
        if (alphav < fabs(v->data[v->row * i + 1 ]))
        {
            alphav = fabs(v->data[v->row * i + 1]);
        }
        if (alphav < fabs(v->data[v->row * i + 2]))
        {
            alphav = fabs(v->data[v->row * i + 2]);
        }
    }
    alphav = alphav / Vmax;

    //alphaa = (max(abs([a(1, :) a(2, :) a(3, :)])) / Amax) ^ (1 / 2);
    double alphaa;
    alphaa = 0.00;
    for (i = 0; i < a->column; i++)
    {
        if (alphaa < fabs(a->data[a->row * i]))
        {
            alphaa = fabs(a->data[a->row * i]);
        }
        if (alphaa < fabs(a->data[a->row * i + 1]))
        {
            alphaa = fabs(a->data[a->row * i + 1]);
        }
        if (alphaa < fabs(a->data[a->row * i + 2]))
        {
            alphaa = fabs(a->data[a->row * i + 2]);
        }
    }
    alphaa = sqrt(alphaa / Amax);

    //alphaj = (max(abs([jj(1, :) jj(2, :) jj(3, :)])) / Jmax) ^ (1 / 3);
    double alphaj;
    alphaj = 0.00;
    for (i = 0; i < jj->column; i++)
    {
        if (alphaj < fabs(jj->data[jj->row * i]))
        {
            alphaj = fabs(jj->data[jj->row * i]);
        }
        if (alphaj < fabs(jj->data[jj->row * i + 1]))
        {
            alphaj = fabs(jj->data[jj->row * i + 1]);
        }
        if (alphaj < fabs(jj->data[jj->row * i + 2]))
        {
            alphaj = fabs(jj->data[jj->row * i + 2]);
        }
    }
    alphaj = pow(alphaj / Jmax, 1.0 / 3);

    //alpharv_A = max(abs(v(4, :))) / RVmax_A;
    double alpharv_A;
    alpharv_A = 0.00;
    for (i = 0; i < v->column; i++)
    {
        if (alpharv_A < fabs(v->data[v->row * i + 3]))
        {
            alpharv_A = fabs(v->data[v->row * i + 3]);
        }

    }
    alpharv_A = alpharv_A / RVmax_A;
    //alpharv_C = max(abs(v(5, :))) / RVmax_C;
    double alpharv_C;
    alpharv_C = 0.00;
    for (i = 0; i < v->column; i++)
    {
        if (alpharv_C < fabs(v->data[v->row * i + 4]))
        {
            alpharv_C = fabs(v->data[v->row * i + 4]);
        }

    }
    alpharv_C = alpharv_C / RVmax_C;

    //alphara = (max(abs([a(4, :) a(5, :)])) / RAmax) ^ (1 / 2);
    double alphara;
    alphara = 0.00;
    for (i = 0; i < a->column; i++)
    {
        if (alphara < fabs(a->data[a->row * i + 3]))
        {
            alphara = fabs(a->data[a->row * i + 3]);
        }
        if (alphara < fabs(a->data[a->row * i + 4]))
        {
            alphara = fabs(a->data[a->row * i + 4]);
        }
    }
    alphara = pow(alphara / RAmax, 1.0 / 2);

    //alpharj = (max(abs([jj(4, :) jj(5, :)])) / RJmax) ^ (1 / 3);
    double alpharj;
    alpharj = 0.00;
    for (i = 0; i < jj->column; i++)
    {
        if (alpharj < fabs(jj->data[jj->row * i + 3]))
        {
            alpharj = fabs(jj->data[jj->row * i + 3]);
        }
        if (alpharj < fabs(jj->data[jj->row * i + 4]))
        {
            alpharj = fabs(jj->data[jj->row * i + 4]);
        }
    }
    alpharj = pow(alpharj / RJmax, 1.0 / 3);

    printf("globaltrans Func  First Values: \n");
    printf("alphav = %lf\n", alphav);
    printf("alphaa = %lf\n", alphaa);
    printf("alphaj = %lf\n", alphaj);
    printf("alpharv_A = %lf\n", alpharv_A);
    printf("alpharv_C = %lf\n", alpharv_C);
    printf("alphara = %lf\n", alphara);
    printf("alpharj = %lf\n", alpharj);

    M_free(jj);
    M_free(a);
    M_free(v);

    //M_free(tknot1);
    M_free(talpha);
    M_free(ttemp);
    //M_free(tknot);

    M_free(nalpha);
    M_free(nbeta);
    //M_free(ncoef);

    free(t.data);

    free(sp1->coefs);
    free(sp1->knots);
    free(sp1);

    free(sp2->coefs);
    free(sp2->knots);
    free(sp2);

    free(sp3->coefs);
    free(sp3->knots);
    free(sp3);

    free(sp4->coefs);
    free(sp4->knots);
    free(sp4);

    free(sp5->coefs);
    free(sp5->knots);
    free(sp5);

    free(sp6->coefs);
    free(sp6->knots);
    free(sp6);

    M_free(dsp->coefs);
    M_free(dsp->knots);
    free(dsp);

    M_free(ddsp->coefs);
    M_free(ddsp->knots);
    free(ddsp);

    M_free(d3sp->coefs);
    M_free(d3sp->knots);
    free(d3sp);


    //while max([alphav alphaa alphaj alpharv_A alpharv_C alphara alpharj]) > 1
    if ((alphav > 1) || (alphaa > 1) || (alphaj > 1) || (alpharv_A > 1) || (alpharv_C > 1) || (alphara > 1) || (alpharj > 1))
    {
        stusp1* nsp7, * nsp8, * nsp9, * nsp10, *sp7, *sp8, *sp9, *sp10;

        nsp7  = (stusp1*)malloc(sizeof(stusp1));
        nsp8  = (stusp1*)malloc(sizeof(stusp1));
        nsp9  = (stusp1*)malloc(sizeof(stusp1));
        nsp10 = (stusp1*)malloc(sizeof(stusp1));

        sp7  = (stusp1*)malloc(sizeof(stusp1));
        sp8  = (stusp1*)malloc(sizeof(stusp1));
        sp9  = (stusp1*)malloc(sizeof(stusp1));
        sp10 = (stusp1*)malloc(sizeof(stusp1));

        nsp7->number        = nsp->number;
        nsp7->coefs         = (Matrix*)malloc(sizeof(Matrix));
        nsp7->coefs->column = nsp->coefs->column;
        nsp7->coefs->row    = nsp->coefs->row;
        nsp7->coefs->data   = (double*)malloc((nsp7->coefs->column * nsp7->coefs->row) * sizeof(double));
        for (i = 0; i < (nsp7->coefs->column * nsp7->coefs->row); i++)
        {
            nsp7->coefs->data[i] = nsp->coefs->data[i];
        }
        nsp7->knots         = (Matrix*)malloc(sizeof(Matrix));
        nsp7->knots->column = nsp->knots->column;
        nsp7->knots->row    = nsp->knots->row;
        nsp7->knots->data   = (double*)malloc((nsp7->knots->column * nsp7->knots->row) * sizeof(double));
        for (i = 0; i < (nsp7->knots->column * nsp7->knots->row); i++)
        {
            nsp7->knots->data[i] = nsp->knots->data[i];
        }

        nsp8->number        = nsp->number;
        nsp8->coefs         = (Matrix*)malloc(sizeof(Matrix));
        nsp8->coefs->column = nsp->coefs->column;
        nsp8->coefs->row    = nsp->coefs->row;
        nsp8->coefs->data   = (double*)malloc((nsp8->coefs->column * nsp8->coefs->row) * sizeof(double));
        for (i = 0; i < (nsp8->coefs->column * nsp8->coefs->row); i++)
        {
            nsp8->coefs->data[i] = nsp->coefs->data[i];
        }
        nsp8->knots         = (Matrix*)malloc(sizeof(Matrix));
        nsp8->knots->column = nsp->knots->column;
        nsp8->knots->row    = nsp->knots->row;
        nsp8->knots->data   = (double*)malloc((nsp8->knots->column * nsp8->knots->row) * sizeof(double));
        for (i = 0; i < (nsp8->knots->column * nsp8->knots->row); i++)
        {
            nsp8->knots->data[i] = nsp->knots->data[i];
        }

        nsp9->number        = nsp->number;
        nsp9->coefs         = (Matrix*)malloc(sizeof(Matrix));
        nsp9->coefs->column = nsp->coefs->column;
        nsp9->coefs->row    = nsp->coefs->row;
        nsp9->coefs->data   = (double*)malloc((nsp9->coefs->column * nsp9->coefs->row) * sizeof(double));
        for (i = 0; i < (nsp9->coefs->column * nsp9->coefs->row); i++)
        {
            nsp9->coefs->data[i] = nsp->coefs->data[i];
        }
        nsp9->knots         = (Matrix*)malloc(sizeof(Matrix));
        nsp9->knots->column = nsp->knots->column;
        nsp9->knots->row    = nsp->knots->row;
        nsp9->knots->data   = (double*)malloc((nsp9->knots->column * nsp9->knots->row) * sizeof(double));
        for (i = 0; i < (nsp9->knots->column * nsp9->knots->row); i++)
        {
            nsp9->knots->data[i] = nsp->knots->data[i];
        }

        nsp10->number        = nsp->number;
        nsp10->coefs         = (Matrix*)malloc(sizeof(Matrix));
        nsp10->coefs->column = nsp->coefs->column;
        nsp10->coefs->row    = nsp->coefs->row;
        nsp10->coefs->data   = (double*)malloc((nsp10->coefs->column * nsp10->coefs->row) * sizeof(double));
        for (i = 0; i < (nsp10->coefs->column * nsp10->coefs->row); i++)
        {
            nsp10->coefs->data[i] = nsp->coefs->data[i];
        }
        nsp10->knots         = (Matrix*)malloc(sizeof(Matrix));
        nsp10->knots->column = nsp->knots->column;
        nsp10->knots->row    = nsp->knots->row;
        nsp10->knots->data   = (double*)malloc((nsp10->knots->column * nsp10->knots->row) * sizeof(double));
        for (i = 0; i < (nsp10->knots->column * nsp10->knots->row); i++)
        {
            nsp10->knots->data[i] = nsp->knots->data[i];
        }

        sp7->number        = sp->number;
        sp7->coefs         = (Matrix*)malloc(sizeof(Matrix));
        sp7->coefs->column = sp->coefs->column;
        sp7->coefs->row    = sp->coefs->row;
        sp7->coefs->data   = (double*)malloc((sp7->coefs->column * sp7->coefs->row) * sizeof(double));
        for (i = 0; i < (sp7->coefs->column * sp7->coefs->row); i++)
        {
            sp7->coefs->data[i] = sp->coefs->data[i];
        }
        sp7->knots         = (Matrix*)malloc(sizeof(Matrix));
        sp7->knots->column = sp->knots->column;
        sp7->knots->row    = sp->knots->row;
        sp7->knots->data   = (double*)malloc((sp7->knots->column * sp7->knots->row) * sizeof(double));
        for (i = 0; i < (sp7->knots->column * sp7->knots->row); i++)
        {
            sp7->knots->data[i] = sp->knots->data[i];
        }



        //[~, alpha1, ~] = bangcoef(nsp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);
        Matrix* ntknot1, * nalpha1, * nttemp, * ntknot;

        ntknot1 = (Matrix*)malloc(sizeof(Matrix));
        nalpha1 = (Matrix*)malloc(sizeof(Matrix));
        nttemp  = (Matrix*)malloc(sizeof(Matrix));

        //talpha = 2 * ones(1, n - 3);
        nalpha1  = Matrix_gen1(1, nsp7->number - 3, 1 * 2);
        //ttemp  = zeros(1, n - 3);
        nttemp   = Matrix_gen1(1, nsp7->number - 3, 0);

        bangcoef(ntknot1, nalpha1, nttemp, nsp7, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);

        //[tknot, alpha, beta] = multal(min(alpha1), nsp);
        double nalpha;
        Matrix * nnalpha, * nnbeta, * nncoef;
        //alpha = min(talpha);

        i = 0;
        nalpha = nalpha1->data[0];
        for (i = 0; i < nalpha1->column; i++)
        {
            if (nalpha > nalpha1->data[i])
            {
                nalpha = nalpha1->data[i];
            }
        }

        nnalpha = Matrix_gen1(1, nsp8->number - 3, 1 * nalpha);
        nnbeta  = Matrix_gen1(1, nsp8->number - 3, 0);
        ntknot  = Matrix_gen1(1, nsp8->number + 4, 0);

        multal(ntknot, nnalpha, nnbeta, nalpha, nsp8);

        //ncoef = bangc(nsp, tknot, alpha, beta);
        nncoef  = Matrix_gen1(5, nsp9->number, 0);
        bangc(nncoef, *nsp9, *ntknot, *nnalpha, *nnbeta);
        
        //nsp = spmak2(tknot, ncoef);
        nsp->number = nncoef->column;
        nsp->knots  = ntknot;
        nsp->coefs  = nncoef;

        //sp = nsp;
        //gap = sp.knots(end) / 10000;
        sp = nsp;
        gap = sp->knots->data[sp->knots->column-1] / 10000;

        sp8->number        = nsp->number;
        sp8->coefs         = (Matrix*)malloc(sizeof(Matrix));
        sp8->coefs->column = nsp->coefs->column;
        sp8->coefs->row    = nsp->coefs->row;
        sp8->coefs->data   = (double*)malloc((sp8->coefs->column * sp8->coefs->row) * sizeof(double));
        for (i = 0; i < (sp8->coefs->column * sp8->coefs->row); i++)
        {
            sp8->coefs->data[i] = nsp->coefs->data[i];
        }
        sp8->knots         = (Matrix*)malloc(sizeof(Matrix));
        sp8->knots->column = nsp->knots->column;
        sp8->knots->row    = nsp->knots->row;
        sp8->knots->data   = (double*)malloc((sp8->knots->column * sp8->knots->row) * sizeof(double));
        for (i = 0; i < (sp8->knots->column * sp8->knots->row); i++)
        {
            sp8->knots->data[i] = nsp->knots->data[i];
        }

        sp9->number        = nsp->number;
        sp9->coefs         = (Matrix*)malloc(sizeof(Matrix));
        sp9->coefs->column = nsp->coefs->column;
        sp9->coefs->row    = nsp->coefs->row;
        sp9->coefs->data   = (double*)malloc((sp9->coefs->column * sp9->coefs->row) * sizeof(double));
        for (i = 0; i < (sp9->coefs->column * sp9->coefs->row); i++)
        {
            sp9->coefs->data[i] = nsp->coefs->data[i];
        }
        sp9->knots         = (Matrix*)malloc(sizeof(Matrix));
        sp9->knots->column = nsp->knots->column;
        sp9->knots->row    = nsp->knots->row;
        sp9->knots->data   = (double*)malloc((sp9->knots->column * sp9->knots->row) * sizeof(double));
        for (i = 0; i < (sp9->knots->column * sp9->knots->row); i++)
        {
            sp9->knots->data[i] = nsp->knots->data[i];
        }

        sp10->number        = nsp->number;
        sp10->coefs         = (Matrix*)malloc(sizeof(Matrix));
        sp10->coefs->column = nsp->coefs->column;
        sp10->coefs->row    = nsp->coefs->row;
        sp10->coefs->data   = (double*)malloc((sp10->coefs->column * sp10->coefs->row) * sizeof(double));
        for (i = 0; i < (sp10->coefs->column * sp10->coefs->row); i++)
        {
            sp10->coefs->data[i] = nsp->coefs->data[i];
        }
        sp10->knots         = (Matrix*)malloc(sizeof(Matrix));
        sp10->knots->column = nsp->knots->column;
        sp10->knots->row    = nsp->knots->row;
        sp10->knots->data   = (double*)malloc((sp10->knots->column * sp10->knots->row) * sizeof(double));
        for (i = 0; i < (sp10->knots->column * sp10->knots->row); i++)
        {
            sp10->knots->data[i] = nsp->knots->data[i];
        }

        //t = 0:gap:sp.knots(end);
        int i;
        i = 0;
        Matrix  t;
        t.data = (double*)malloc(10001 * sizeof(double));
        t.row = 1;
        t.column = 10001;
        for (i = 0; i <= 10000; i++)
        {
            t.data[i] = gap * i; 
        }

        //stusp1* dsp, * ddsp, * d3sp;
        //Matrix* v, * a, * jj;
        dsp  = (stusp1*)malloc(sizeof(stusp1));
        ddsp = (stusp1*)malloc(sizeof(stusp1));
        d3sp = (stusp1*)malloc(sizeof(stusp1));

        //% 速度
        //dsp = fnder2(sp, 1);
        //v = fnval2(dsp, t);
        fnder2(dsp, sp8, 1);
        v = Matrix_gen1(sp8->coefs->row, t.column, 0.00);
        fnval2(v, dsp, t.data, t.column);

         //% 加速度
         //ddsp = fnder2(sp, 2);
         //a = fnval2(ddsp, t);
        fnder2(ddsp, sp9, 2);
        a = Matrix_gen1(ddsp->coefs->row, t.column, 0.00);
        fnval2(a, ddsp, t.data, t.column);

         //% 加加速
         //d3sp = fnder2(sp, 3);
         //jj = fnval2(d3sp, t);
        fnder2(d3sp, sp10, 3);
        jj = Matrix_gen1(d3sp->coefs->row, t.column, 0.00);
        fnval2(jj, d3sp, t.data, t.column);

 
        //alphav = max(abs([v(1, :) v(2, :) v(3, :)])) / Vmax;
        //double alphav;
        alphav = 0.00;
        for (i = 0; i < v->column; i++)
        {
            if (alphav < fabs(v->data[v->row * i]))
            {
                alphav = fabs(v->data[v->row * i]);
            }
            if (alphav < fabs(v->data[v->row * i + 1]))
            {
                alphav = fabs(v->data[v->row * i + 1]);
            }
            if (alphav < fabs(v->data[v->row * i + 2]))
            {
                alphav = fabs(v->data[v->row * i + 2]);
            }
        }
        alphav = alphav / Vmax;

        //alphaa = (max(abs([a(1, :) a(2, :) a(3, :)])) / Amax) ^ (1 / 2);
        //double alphaa;
        alphaa = 0.00;
        for (i = 0; i < a->column; i++)
        {
            if (alphaa < fabs(a->data[a->row * i]))
            {
                alphaa = fabs(a->data[a->row * i]);
            }
            if (alphaa < fabs(a->data[a->row * i + 1]))
            {
                alphaa = fabs(a->data[a->row * i + 1]);
            }
            if (alphaa < fabs(a->data[a->row * i + 2]))
            {
                alphaa = fabs(a->data[a->row * i + 2]);
            }
        }
        alphaa = sqrt(alphaa / Amax);

        //alphaj = (max(abs([jj(1, :) jj(2, :) jj(3, :)])) / Jmax) ^ (1 / 3);
        //double alphaj;
        alphaj = 0.00;
        for (i = 0; i < jj->column; i++)
        {
            if (alphaj < fabs(jj->data[jj->row * i]))
            {
                alphaj = fabs(jj->data[jj->row * i]);
            }
            if (alphaj < fabs(jj->data[jj->row * i + 1]))
            {
                alphaj = fabs(jj->data[jj->row * i + 1]);
            }
            if (alphaj < fabs(jj->data[jj->row * i + 2]))
            {
                alphaj = fabs(jj->data[jj->row * i + 2]);
            }
        }
        alphaj = pow(alphaj / Jmax, 1.0 / 3);

        //alpharv_A = max(abs(v(4, :))) / RVmax_A;
        //double alpharv_A;
        alpharv_A = 0.00;
        for (i = 0; i < v->column; i++)
        {
            if (alpharv_A < fabs(v->data[v->row * i + 3]))
            {
                alpharv_A = fabs(v->data[v->row * i + 3]);
            }

        }
        alpharv_A = alpharv_A / RVmax_A;
        //alpharv_C = max(abs(v(5, :))) / RVmax_C;
        //double alpharv_C;
        alpharv_C = 0.00;
        for (i = 0; i < v->column; i++)
        {
            if (alpharv_C < fabs(v->data[v->row * i + 4]))
            {
                alpharv_C = fabs(v->data[v->row * i + 4]);
            }

        }
        alpharv_C = alpharv_C / RVmax_C;

        //alphara = (max(abs([a(4, :) a(5, :)])) / RAmax) ^ (1 / 2);
        //double alphara;
        alphara = 0.00;
        for (i = 0; i < a->column; i++)
        {
            if (alphara < fabs(a->data[a->row * i + 3]))
            {
                alphara = fabs(a->data[a->row * i + 3]);
            }
            if (alphara < fabs(a->data[a->row * i + 4]))
            {
                alphara = fabs(a->data[a->row * i + 4]);
            }
        }
        alphara = pow(alphara / RAmax, 1.0 / 2);

        //alpharj = (max(abs([jj(4, :) jj(5, :)])) / RJmax) ^ (1 / 3);
        //double alpharj;
        alpharj = 0.00;
        for (i = 0; i < jj->column; i++)
        {
            if (alpharj < fabs(jj->data[jj->row * i + 3]))
            {
                alpharj = fabs(jj->data[jj->row * i + 3]);
            }
            if (alpharj < fabs(jj->data[jj->row * i + 4]))
            {
                alpharj = fabs(jj->data[jj->row * i + 4]);
            }
        }
        alpharj = pow(alpharj / RJmax, 1.0 / 3);


        printf("globaltrans Func :: while max([alphav alphaa alphaj alpharv_A alpharv_C alphara alpharj]) >1, Second Values:\n");
        printf("alphav = %lf\n", alphav);
        printf("alphaa = %lf\n", alphaa);
        printf("alphaj = %lf\n", alphaj);
        printf("alpharv_A = %lf\n", alpharv_A);
        printf("alpharv_C = %lf\n", alpharv_C);
        printf("alphara = %lf\n", alphara);
        printf("alpharj = %lf\n", alpharj);

        M_free(jj);
        M_free(a);
        M_free(v);
        free(t.data);


        M_free(nnalpha);
        M_free(nnbeta);
        M_free(nncoef);

        //M_free(ntknot);
        M_free(ntknot);
        //M_free(ntknot1);
        M_free(nalpha1);
        M_free(nttemp);
        M_free(tknot);

        //M_free(nbeta);
        M_free(ncoef);

        free(sp7->coefs);
        free(sp7->knots);
        free(sp7);

        free(sp8->coefs);
        free(sp8->knots);
        free(sp8);

        free(sp9->coefs);
        free(sp9->knots);
        free(sp9);

        free(sp10->coefs);
        free(sp10->knots);
        free(sp10);

        free(nsp7->coefs);
        free(nsp7->knots);
        free(nsp7);

        free(nsp8->coefs);
        free(nsp8->knots);
        free(nsp8);

        free(nsp9->coefs);
        free(nsp9->knots);
        free(nsp9);

        free(nsp10->coefs);
        free(nsp10->knots);
        free(nsp10);


        M_free(dsp->coefs);
        M_free(dsp->knots);
        free(dsp);

        M_free(ddsp->coefs);
        M_free(ddsp->knots);
        free(ddsp);

        M_free(d3sp->coefs);
        M_free(d3sp->knots);
        free(d3sp);

    }

}

