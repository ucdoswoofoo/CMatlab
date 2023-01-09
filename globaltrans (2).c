#include "pch.h"
#include "globaltrans.h"
#include "bangcoef.h"
#include "multal.h"
#include "bangc.h"
#include "spmak2.h"
#include "fnder2.h"
#include "fnval2.h"

/*
% 整体调整曲线节点，使得曲线形状不变，时间延长或缩短, 曲线满足运动约束
% 输入：1 sp 样条结构体 2 Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax double型1 * 1
% 输出：nsp 样条结构体
function nsp = globaltrans(sp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax)
*/


void globaltrans(stusp* nsp, const stusp* sp, const double  Vmax, const double  Amax, const double  Jmax, const double RVmax_A, const double RVmax_C, const double  RAmax, const double  RJmax)
{
    /*
    [~, alpha1, ~] = bangcoef(sp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);
    [tknot, alpha, beta] = multal(min(alpha1), sp);
    ncoef = bangc(sp, tknot, alpha, beta);
    nsp = spmak2(tknot, ncoef);
    */
    /*
    gap = sp.knots(end) / 10000;
    t = 0:gap:sp.knots(end);
    //% 速度
        dsp = fnder2(sp, 1);
    v = fnval2(dsp, t);
    //% 加速度
        ddsp = fnder2(sp, 2);
    a = fnval2(ddsp, t);
    //% 加加速
        d3sp = fnder2(sp, 3);
    jj = fnval2(d3sp, t);
    alphav = max(abs([v(1, :) v(2, :) v(3, :)])) / Vmax;
    alphaa = (max(abs([a(1, :) a(2, :) a(3, :)])) / Amax) ^ (1 / 2);
    alphaj = (max(abs([jj(1, :) jj(2, :) jj(3, :)])) / Jmax) ^ (1 / 3);
    alpharv_A = max(abs(v(4, :))) / RVmax_A;
    alpharv_C = max(abs(v(5, :))) / RVmax_C;
    alphara = (max(abs([a(4, :) a(5, :)])) / RAmax) ^ (1 / 2);
    alpharj = (max(abs([jj(4, :) jj(5, :)])) / RJmax) ^ (1 / 3);
    */
    
    //[~, alpha1, ~] = bangcoef(sp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);
    CellMat tknot, talpha, ttemp;
    bangcoef(tknot, talpha, ttemp,  sp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);

    //[tknot, alpha, beta] = multal(min(alpha1), sp);
    //ncoef = bangc(sp, tknot, alpha, beta);
    //nsp = spmak2(tknot, ncoef);

    CellMat nalpha, nbeta;
    multal(tknot, nalpha, nbeta, talpha.CellData[0][0], sp);

    //ncoef = bangc(sp, tknot, alpha, beta);

    CellMat ncoef;

    bangc( ncoef, sp, tknot, talpha, nbeta);

    //nsp = spmak2(tknot, ncoef);
    spmak2(tknot, ncoef, nsp);


    //gap = sp.knots(end) / 10000;

    double gap;
    gap = sp->knots.CellData[0][sp->knots.iCol] / 10000;

    //t = 0:gap:sp.knots(end);

    int i;
    i = 0;
    CellMat  t;
    //赋值未决
    for (i=0; i<= sp->knots.CellData[0][sp->knots.iCol] * 10000; i++)
    {
        t.CellData[0][i] = sp->knots.CellData[0][i];

    }


    //% 速度
    //dsp = fnder2(sp, 1);

    stusp* dsp,* ddsp, * d3sp;
    fnder2(dsp, sp, 1);

    //v = fnval2(dsp, t);
    CellMat v, a;
    fnval2(v, dsp, t);
    //% 加速度
    //ddsp = fnder2(sp, 2);
    fnder2(ddsp, sp, 2);

    //a = fnval2(ddsp, t);
    fnval2(a, ddsp, t);

    //% 加加速
    //d3sp = fnder2(sp, 3);
    fnder2(d3sp, sp, 3);

    //jj = fnval2(d3sp, t);
    CellMat jj;
    fnval2(jj, d3sp, t);

    //alphav = max(abs([v(1, :) v(2, :) v(3, :)])) / Vmax;

    double alphav;

    alphav = 0.00;
    for (i=0; i< v.iCol; i++)
    {
        if (alphav < abs(v.CellData[0][i]))
        {
            alphav = abs(v.CellData[0][i]);
        }
        if (alphav < abs(v.CellData[1][i]))
        {
            alphav = abs(v.CellData[1][i]);
        }
        if (alphav < abs(v.CellData[2][i]))
        {
            alphav = abs(v.CellData[2][i]);
        }
    }
    alphav = alphav / Vmax;

    //alphaa = (max(abs([a(1, :) a(2, :) a(3, :)])) / Amax) ^ (1 / 2);

    double alphaa;
    alphaa = 0.00;
    for (i = 0; i < a.iCol; i++)
    {
        if (alphaa < abs(a.CellData[0][i]))
        {
            alphaa = abs(a.CellData[0][i]);
        }
        if (alphaa < abs(a.CellData[1][i]))
        {
            alphaa = abs(a.CellData[1][i]);
        }
        if (alphaa < abs(a.CellData[2][i]))
        {
            alphaa = abs(a.CellData[2][i]);
        }
    }

    alphaa = sqrt(alphaa / Amax);



    //alphaj = (max(abs([jj(1, :) jj(2, :) jj(3, :)])) / Jmax) ^ (1 / 3);

    double alphaj;
    alphaj = 0.00;
    for (i = 0; i < jj.iCol; i++)
    {
        if (alphaj < abs(jj.CellData[0][i]))
        {
            alphaj = abs(jj.CellData[0][i]);
        }
        if (alphaj < abs(jj.CellData[1][i]))
        {
            alphaj = abs(jj.CellData[1][i]);
        }
        if (alphaj < abs(jj.CellData[2][i]))
        {
            alphaj = abs(jj.CellData[2][i]);
        }
    }

    alphaj = pow(alphaj / Jmax, 1/3);

    //alpharv_A = max(abs(v(4, :))) / RVmax_A;

    double alpharv_A;
    alpharv_A = 0.00;
    for (i = 0; i < v.iCol; i++)
    {
        if (alpharv_A < abs(v.CellData[3][i]))
        {
            alpharv_A = abs(v.CellData[3][i]);
        }

    }
    alpharv_A = alpharv_A / RVmax_A;

    //alpharv_C = max(abs(v(5, :))) / RVmax_C;

    double alpharv_C;
    alpharv_C = 0.00;
    for (i = 0; i < v.iCol; i++)
    {
        if (alpharv_C < abs(v.CellData[4][i]))
        {
            alpharv_C = abs(v.CellData[4][i]);
        }

    }
    alpharv_C = alpharv_C / RVmax_C;

    //alphara = (max(abs([a(4, :) a(5, :)])) / RAmax) ^ (1 / 2);

    double alphara;
    alphara = 0.00;
    for (i = 0; i < a.iCol; i++)
    {
        if (alphara < abs(a.CellData[3][i]))
        {
            alphara = abs(a.CellData[3][i]);
        }
        if (alphara < abs(a.CellData[4][i]))
        {
            alphara = abs(a.CellData[4][i]);
        }
    }

    alphara = pow(alphara / RAmax, 1 / 2);

    //alpharj = (max(abs([jj(4, :) jj(5, :)])) / RJmax) ^ (1 / 3);
    double alpharj;
    alpharj = 0.00;
    for (i = 0; i < jj.iCol; i++)
    {
        if (alpharj < abs(jj.CellData[3][i]))
        {
            alpharj = abs(jj.CellData[3][i]);
        }
        if (alpharj < abs(jj.CellData[4][i]))
        {
            alpharj = abs(jj.CellData[4][i]);
        }
    }

    alpharj = pow(alpharj / RJmax, 1 / 3);

    //while max([alphav alphaa alphaj alpharv_A alpharv_C alphara alpharj]) > 1
    if ( (alphav>1) || (alphaa>1) || (alphaj>1) || (alpharv_A>1) || (alpharv_C>1) || (alphara >1) || (alpharj>1) )
    {
            //[~, alpha1, ~] = bangcoef(nsp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);

            CellMat tknot, talpha, ttemp, alpha1;
            bangcoef(tknot, alpha1, ttemp, nsp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);


            //[tknot, alpha, beta] = multal(min(alpha1), nsp);

            CellMat alpha,beta;
            multal(tknot, alpha, beta, alpha1.CellData[0][0], nsp);


            //ncoef = bangc(nsp, tknot, alpha, beta);
            //nsp = spmak2(tknot, ncoef);

            CellMat ncoef;

            bangc(ncoef, nsp, tknot, alpha, beta);
            spmak2(tknot, ncoef, nsp);

            //sp = nsp;
            //gap = sp.knots(end) / 10000;
            sp = nsp;
            gap = sp->knots.CellData[0][sp->knots.iCol] / 10000;


            //t = 0:gap:sp.knots(end);

            int i;
            i = 0;
            CellMat  t;
            //赋值未决
            for (i = 0; i <= sp->knots.CellData[0][sp->knots.iCol] * 10000; i++)
            {
                t.CellData[0][i] = sp->knots.CellData[0][i];
            }

            //% 速度
            //dsp = fnder2(sp, 1);
            //v = fnval2(dsp, t);

            stusp* dsp, * ddsp, * d3sp;
            fnder2(dsp, sp, 1);

            //v = fnval2(dsp, t);
            CellMat v, a;
            fnval2(v, dsp, t);

            //% 加速度
            //ddsp = fnder2(sp, 2);
            //a = fnval2(ddsp, t);

            fnder2(ddsp, sp, 2);
            //v = fnval2(dsp, t);
            fnval2(a, ddsp, t);

            //% 加加速
            //d3sp = fnder2(sp, 3);
            //jj = fnval2(d3sp, t);

            fnder2(d3sp, sp, 3);
            //v = fnval2(dsp, t);
            fnval2(jj, d3sp, t);

            //alphav = max(abs([v(1, :) v(2, :) v(3, :)])) / Vmax;
            double alphav;

            alphav = 0.00;
            for (i = 0; i < v.iCol; i++)
            {
                if (alphav < abs(v.CellData[0][i]))
                {
                    alphav = abs(v.CellData[0][i]);
                }
                if (alphav < abs(v.CellData[1][i]))
                {
                    alphav = abs(v.CellData[1][i]);
                }
                if (alphav < abs(v.CellData[2][i]))
                {
                    alphav = abs(v.CellData[2][i]);
                }
            }
            alphav = alphav / Vmax;

            //alphaa = (max(abs([a(1, :) a(2, :) a(3, :)])) / Amax) ^ (1 / 2);
            double alphaa;

            alphaa = 0.00;
            for (i = 0; i < a.iCol; i++)
            {
                if (alphaa < abs(a.CellData[0][i]))
                {
                    alphaa = abs(a.CellData[0][i]);
                }
                if (alphaa < abs(a.CellData[1][i]))
                {
                    alphaa = abs(a.CellData[1][i]);
                }
                if (alphaa < abs(a.CellData[2][i]))
                {
                    alphaa = abs(a.CellData[2][i]);
                }
            }
            alphaa = pow(alphaa / Amax, 1/2);

            //alphaj = (max(abs([jj(1, :) jj(2, :) jj(3, :)])) / Jmax) ^ (1 / 3);

            double alphaj;

            alphaj = 0.00;
            for (i = 0; i < jj.iCol; i++)
            {
                if (alphaj < abs(jj.CellData[0][i]))
                {
                    alphaj = abs(jj.CellData[0][i]);
                }
                if (alphaj < abs(jj.CellData[1][i]))
                {
                    alphaj = abs(jj.CellData[1][i]);
                }
                if (alphaj < abs(jj.CellData[2][i]))
                {
                    alphaj = abs(jj.CellData[2][i]);
                }
            }
            alphaj = pow(alphaj / Jmax, 1 / 3);

            //alpharv_A = max(abs(v(4, :))) / RVmax_A;
            double alpharv_A;

            alpharv_A = 0.00;
            for (i = 0; i < v.iCol; i++)
            {
                if (alpharv_A < abs(v.CellData[3][i]))
                {
                    alpharv_A = abs(v.CellData[3][i]);
                }
            }
            alpharv_A = alpharv_A / RVmax_A;


            //alpharv_C = max(abs(v(5, :))) / RVmax_C;
            double alpharv_C;

            alpharv_C = 0.00;
            for (i = 0; i < v.iCol; i++)
            {
                if (alpharv_C < abs(v.CellData[4][i]))
                {
                    alpharv_C = abs(v.CellData[4][i]);
                }
            }
            alpharv_C = alpharv_C / RVmax_C;


            //alphara = (max(abs([a(4, :) a(5, :)])) / RAmax) ^ (1 / 2);
            double alphara;

            alphara = 0.00;
            for (i = 0; i < a.iCol; i++)
            {
                if (alphara < abs(a.CellData[3][i]))
                {
                    alphara = abs(a.CellData[3][i]);
                }
                if (alphara < abs(a.CellData[4][i]))
                {
                    alphara = abs(a.CellData[4][i]);
                }
            }
            alphara = pow(alphara / RAmax, 1/2);


            //alpharj = (max(abs([jj(4, :) jj(5, :)])) / RJmax) ^ (1 / 3);
            double alpharj;

            alpharj = 0.00;
            for (i = 0; i < jj.iCol; i++)
            {
                if (alpharj < abs(jj.CellData[3][i]))
                {
                    alpharj = abs(jj.CellData[3][i]);
                }
                if (alpharj < abs(jj.CellData[4][i]))
                {
                    alpharj = abs(jj.CellData[4][i]);
                }
            }
            alpharj = pow(alpharj / RJmax, 1 / 3);


    }
    //end




}


