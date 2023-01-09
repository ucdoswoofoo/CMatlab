#include "data3curve.h"
#include "fnval2.h"
#include "fnder2.h"

/*
% 点到插值曲线的距离

function[d, dmj, fpt] = data3curve(segwrefq, sp)
r = length(segwrefq(1, :));

fpt = zeros(r, 1);
gap = 0.01;
sampfpt = 0:gap:1;
sampfpt = (sampfpt - sp.knots(1)).*(sp.knots(end) - sp.knots(1));
samp = fnval2(sp, sampfpt); % matlab自带函数，可用自定义函数fnval2替代

*/



void  data3curve(double d, double  dmj, Matrix fpt, Matrix* segwrefq, const stusp1* sp)
{
    //r = length(segwrefq(1, :));
    /*
    int r;
    
    r = 0;
    r = segwrefq->iCol;

    //fpt = zeros(r, 1);
    ZXMatInit(fpt, r, r, 1);

    //gap = 0.01;
    //sampfpt = 0:gap:1;

    double gap;
    int i;
    CellMat sampfpt;

    i = 0;
    gap = 0.01;
    sampfpt.iRow = 1;
    sampfpt.iCol = 100 / 1; // gap;

    for (i=0; i < sampfpt.iCol; i++)
    {
        sampfpt.CellData[0][i] = gap * i;
    }

    //sampfpt = (sampfpt - sp.knots(1)).*(sp.knots(end) - sp.knots(1));
    //samp = fnval2(sp, sampfpt); % matlab自带函数，可用自定义函数fnval2替代
    for (i = 0; i < sampfpt.iCol; i++)
    {
        sampfpt.CellData[0][i] = (sampfpt.CellData[0][i] - sp->knots.CellData[0][0]) * (sp->knots.CellData[0][sp->knots.iCol] - sp->knots.CellData[0][0]);
    }

    CellMat samp, disx, disy, disz, dis;

   // fnval2(samp, sp, sampfpt);

    
    //for i = 1:r
    //    disx = samp(1, :) - segwrefq(1, i);
    //    disy = samp(2, :) - segwrefq(2, i);
    //    disz = samp(3, :) - segwrefq(3, i);
    //    dis = disx. ^ 2 + disy. ^ 2 + disz. ^ 2;
    //    [mm, mi] = min(dis);
    //    fpt(i) = (mi - 1) * gap;
    //end
    

    int j, mi;
    double mm;

    mi = 0;
    mm = 0.00;
    j  = 0;

    for (i=0; i < r; i++)
    {
        for (j=0; j < segwrefq->iCol-1; j++)
        {
            disx.CellData[0][j] = samp.CellData[0][j] - segwrefq->CellData[0][j];
            disy.CellData[0][j] = samp.CellData[1][j] - segwrefq->CellData[1][j];
            disz.CellData[0][j] = samp.CellData[2][j] - segwrefq->CellData[2][j];
        }

        for (j = 0; j < segwrefq->iCol - 1; j++)
        {
            dis.CellData[0][j] = pow(disx.CellData[0][j], 2) + pow(disy.CellData[0][j], 2) + pow(disz.CellData[0][j], 2);
        }

        for (j = 0; j < segwrefq->iCol - 1; j++)
        {
            dis.CellData[0][j] = pow(disx.CellData[0][j], 2) + pow(disy.CellData[0][j], 2) + pow(disz.CellData[0][j], 2);
        }
        for (j = 0; j < segwrefq->iCol - 1; j++)
        {
            if ( 0 == j )
            {
                mm = dis.CellData[0][j];
                mi = 1;
            }

            if ( mm > dis.CellData[0][j] )
            {
                mm = dis.CellData[0][j];
                mi = j + 1;
            }
        }

        fpt.CellData[0][i] = (mi - 1) * gap;

    }
    
    //iter = ones(r, 1);
    //dsp = fnder2(sp, 1); % fnder计算样条曲线的一阶导数
    //ddsp = fnder2(sp, 2); % fnder计算样条曲线的二阶导数
    //step = ones(r, 1) * 2;
    //sig = zeros(r, 1);
    //beta = zeros(r, 1);
    //alpha = zeros(r, 1);
    


    CellMat iter, step, sig, beta, alpha;
    stusp* dsp, * ddsp;

    ZXMatInit(iter, 1, r, 1);

    //fnder2(dsp, sp, 1);
   // fnder2(ddsp, sp, 2);
    ZXMatInit(step, 1, r, 2);
    ZXMatInit(sig, r, r, 1);
    ZXMatInit(beta, r, r, 1);
    ZXMatInit(alpha, r, r, 1);

    double Douiter;

    Douiter = 0.00;

    for (i=0; i< iter.iCol; i++)
    {
        Douiter = Douiter + iter.CellData[0][i];
    }

    while (Douiter > 0)
    {
        CellMat fp, dfp, ddfp, fpda;
        
        //fp = fnval2(sp, fpt);
        //dfp = fnval2(dsp, fpt);
        //ddfp = fnval2(ddsp, fpt);
        
       // fnval2(fp, sp, fpt);
       // fnval2(dfp, dsp, fpt);
       // fnval2(ddfp, ddsp, fpt);

        //fpda = fp - segwrefq;

        for (i=0; i< segwrefq->iRow; i++)
        {
            for (j=0; j< segwrefq->iCol; j++)
            {
                fpda.CellData[i][j] = fp.CellData[i][j] - segwrefq->CellData[i][j];
            }
        }
        //% 这里是一个关于距离的二次函数求极小值，df = 0求得critical点，ddf = df'>0达到极小值
        //f = (fpda(1, :). ^ 2 + fpda(2, :). ^ 2 + fpda(3, :). ^ 2). / 2;

        CellMat f, df, ddf;


        for (i = 0; i < fpda.iCol; i++)
        {
            f.CellData[0][i] = (pow(fpda.CellData[0][i], 2) + pow(fpda.CellData[1][i], 2) + pow(fpda.CellData[2][i], 2)) / 2;
        }

        //df = fpda(1, :).*dfp(1, :) + fpda(2, :).*dfp(2, :) + fpda(3, :).*dfp(3, :);
        //ddf = dfp(1, :). ^ 2 + dfp(2, :). ^ 2 + dfp(3, :). ^ 2 + fpda(1, :).*ddfp(1, :) + fpda(2, :).*ddfp(2, :) + fpda(3, :).*ddfp(3, :);

        for (i = 0; i < fpda.iCol; i++)
        {
            df.CellData[0][i] = fpda.CellData[0][i] * dfp.CellData[0][i] + fpda.CellData[1][i] * dfp.CellData[1][i] + fpda.CellData[2][i] * dfp.CellData[2][i];
        }       
        
        for (i = 0; i < fpda.iCol; i++)
        {
            ddf.CellData[0][i] = pow(dfp.CellData[0][i], 2) + pow(dfp.CellData[1][i], 2) + pow(dfp.CellData[2][i], 2) + fpda.CellData[0][i] * ddfp.CellData[0][i] + fpda.CellData[1][i] * ddfp.CellData[1][i] + fpda.CellData[2][i] * ddfp.CellData[2][i];
        }


        for (i = 1; i <= r; i++)
        {
            if (iter.CellData[0][i] > 0)
            {
                switch ((int)(step.CellData[0][i]))
                {
                    case 2:
                    {
                        
                        //if df(i) >= 10 ^ -4 || df(i) <= -10 ^ -4
                        //    step(i) = 4;
                        //elseif ddf(i) >= 0
                        //    iter(i) = 0;
                        //else
                        //    sig(i) = gap;
                        //step(i) = 3;
                        //end
                        

                        if ( ( df.CellData[0][i] >= pow(10, -4) ) || ( df.CellData[0][i] <= -pow( 10, -4) ))
                        {
                            step.CellData[0][i] = 4;
                        }
                        else if (ddf.CellData[0][i] >= 0)
                        {
                            iter.CellData[0][i] = 0;
                        }
                        else
                        {
                            sig.CellData[0][i] = gap;
                            sig.CellData[i%sig.iRow][i/sig.iRow] = gap;
                            step.CellData[0][i] = 3;
                        }
                    }
                    case 3:
                    {
                        
                        //sig(i) = 2 * sig(i);
                        //if fnval2(sp, fpt(i) + sig(i)) >= f(i)
                        //    step(i) = 3;
                        //else
                        //    fpt(i) = fpt(i) + sig(i);
                        //step(i) = 2;
                        //end
                        

                        sig.CellData[0][i] = 2 * sig.CellData[0][i];

                        CellMat fnval2Tmp, doufnvalsp;

                        doufnvalsp.CellData[0][i] = fpt.CellData[0][i] + sig.CellData[0][i];

                       // fnval2(fnval2Tmp, sp, doufnvalsp);

                        if (fnval2Tmp.CellData[0][i] >= f.CellData[0][i])
                        {
                            step.CellData[0][i] = 3;
                        }
                        else
                        {
                            fpt.CellData[0][i] = fpt.CellData[0][i] + sig.CellData[0][i];
                            step.CellData[0][i] = 2;
                        }
                    }
                    case 4:
                    {
                        
                        //beta(i) = ddf(i);
                        //if beta(i) <= 0
                        //    beta(i) = 1;
                        //end
                        //    alpha(i) = 1;
                        //step(i) = 5;
                        

                        beta.CellData[0][i] = ddf.CellData[0][i];
                        if (beta.CellData[0][i] <= 0)
                        {
                            beta.CellData[0][i] = 1;
                        }

                        alpha.CellData[0][i] = 1;
                        step.CellData[0][i] = 5;
                    }
                    case 5:
                    {
                        
                        //fna = fnval2(sp, fpt(i) - alpha(i) * df(i) / beta(i));
                        //fa = ((fna(1) - segwrefq(1, i)) ^ 2 + (fna(2) - segwrefq(2, i)) ^ 2 + (fna(3) - segwrefq(3, i)) ^ 2) / 2;
                        //fb = f(i) - df(i) ^ 2 * alpha(i) / beta(i) / 4;
                        
                        //if fa <= fb
                        //    step(i) = 6;
                        //else
                        //    alpha(i) = alpha(i) / 2;
                        //    step(i) = 5;
                        //end
                        
                        CellMat fna, fa, fb, TmpMat;

                        TmpMat.CellData[0][i] =   (fpt.CellData[0][i] - alpha.CellData[0][i] * df.CellData[0][i] / beta.CellData[0][i]);
                       // fnval2(fna, sp, TmpMat);

                        fa.CellData[0][i] = (pow((fna.CellData[0][0] - segwrefq->CellData[0][i]), 2) + pow((fna.CellData[0][1] - segwrefq->CellData[1][i]), 2) + pow((fna.CellData[0][2] - segwrefq->CellData[2][i]), 2)) / 2;

                        //fb = f(i) - df(i) ^ 2 * alpha(i) / beta(i) / 4;
                        fb.CellData[0][i] = f.CellData[0][0] - pow(df.CellData[0][i], 2) * alpha.CellData[0][i] / 4;

                        if (fa.CellData[0][i] <= fb.CellData[0][i] )
                        {
                            step.CellData[0][i] = 6;
                        }
                        else
                        {
                            alpha.CellData[0][i] = alpha.CellData[0][i] * 2;
                            step.CellData[0][i] = 5;
                        }
                    }
                    case 6:
                    {
                       //fpt(i) = fpt(i) - alpha(i) * df(i) / beta(i);
                        //ap(i) = 2;


                        fpt.CellData[0][i] = fpt.CellData[0][i] - alpha.CellData[0][i] * df.CellData[0][i] / beta.CellData[0][i];
                        step.CellData[0][i] = 2;

                    }

                }


            }
            
            //if fpt(i) <= sp.knots(1)
            //    fpt(i) = sp.knots(1);
            //iter(i) = 0;
            //elseif fpt(i) >= sp.knots(end)
            //  fpt(i) = sp.knots(end);
            //iter(i) = 0;
            //end
            

            if (fpt.CellData[0][i] <= sp->knots.CellData[0][0])
            {
                fpt.CellData[0][i] = sp->knots.CellData[0][i];
            }
            else if (fpt.CellData[0][i] >= sp->knots.CellData[0][sp->knots.iCol-1])
            {
                fpt.CellData[0][i] = sp->knots.CellData[0][sp->knots.iCol - 1];
                iter.CellData[0][i] = 0;
            }

        }



    }


    //fp = fnval2(sp, fpt);

    CellMat fp, fpda;
    //fnval2(fp, sp, fpt);

    //fpda = fp - segwrefq;

    for (i = 0; i < segwrefq->iRow; i++)
    {
        for (j = 0; j < segwrefq->iCol; j++)
        {
            fpda.CellData[i][j] = fp.CellData[i][j] - segwrefq->CellData[i][j];
        }
    }


    //dis = (fpda(1, :)'.^2+fpda(2,:)'. ^ 2 + fpda(3, :)'.^2).^0.5;


    for (i = 0; i < fpda.iCol; i++)
    {
        dis.CellData[0][i] =  sqrt(pow(fpda.CellData[0][i], 2) + pow(fpda.CellData[1][i], 2) + pow(fpda.CellData[2][i], 2));
    }

        

    //[d, dmj] = max(dis);
    //double d;
    //int dmj;

    for (j = 0; j < dis.iCol - 1; j++)
    {
        if (0 == j)
        {
            d = dis.CellData[0][j];
            dmj = 1;
        }

        if (d < dis.CellData[0][j])
        {
            d = dis.CellData[0][j];
            dmj = j + 1;
        }
    }
    */

}



