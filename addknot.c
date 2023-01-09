#include "addknot.h"
#include "ZXMat.h"
#include "matrix.h"
#include "Binsert.h"
#include "spmak2.h"


void addknot(stusp1* nsp, Matrix* Primsuperknots, const stusp1* sp, Matrix* pkd, double wchorderr, double wdataerr, double eps, Matrix* pkpar, Matrix* dchord, double rd, double reps)
{
    //disp('knot addition');
    //% 添加刀轴点误差超界的极值点
    //superknots1 = zeros(1, length(pkd));
    Matrix  *superknots1, *superknots2, *t, *superknots1Tmp, * ncontr, *Matkk;
    int i, j, m, n, N, kk, iCounter, iMatkkZeroCnt;;
    double gap, * superknots1Data;

    //printf("knot addition\n\r");

    i = 0;
    j = 0;
    m = 0;
    n = 0;
    N = 0;
    iCounter = 0;
    kk = 0;
    gap = 0.00;
    superknots1Data = (double*)malloc(pkd->column * sizeof(double));
    superknots1     = Matrix_gen1(1, pkd->column, 0.00);

    //for i = 1:length(pkd)
    iCounter = 0;
    for (i=0; i < pkd->column; i++)
    {
        //if pkd(i) > max(max(wchorderr, wdataerr) - 0.005, eps)
        if (pkd->data[i] > fmax(fmax(wchorderr, wdataerr) - 0.005, eps))
        {
            //superknots1(i) = pkpar(i);
            superknots1Data[iCounter] = pkpar->data[i];
            superknots1->data[i] = pkpar->data[i];
            iCounter                  = iCounter + 1;
        }
    }

    //superknots1(superknots1 == 0) = [];
    m = 0;
    for (i=0; i< superknots1->column; i++)
    {
        if ( fabs(superknots1->data[i]) > 0.00000001)
        {
            superknots1Data[m] = superknots1->data[i];
            m = m + 1;
        }
    }

    n = 0;
    for (i=0; i< superknots1->column; i++)
    {
        if (fabs(superknots1->data[i]) > 0.00000001)
        {
            superknots1->data[n] = superknots1Data[n];
            n = n + 1;
        }
    }
    superknots1->column = m;
    free(superknots1Data);
    //% 添加刀轴方向误差超界的局部极值点
    //N = 5;
    N = 5;

    //t = 0:sp.knots(end) / N : sp.knots(end);
    gap = sp->knots->data[sp->knots->column-1] / N;
    t = Matrix_gen1(1, N + 1, 0.00);
    t->data[0] = 0.00;
    for (i=0; i<N+1 ; i++)
    {
        t->data[i] = gap * i;
    }

    //superknots2 = zeros(1, 10000);
    superknots2  = Matrix_gen1(1, 10000, 0.00);
    //for i = 2:length(dchord) - 1
    //    if dchord(i) > dchord(i - 1) && dchord(i) > dchord(i + 1) && dchord(i) > max(rd - 0.0005, reps)
    //        superknots2(i) = t(i);
    //   end
    //end
    for (i=2; i< max(dchord->column, dchord->row); i++)
    {
        if ( (dchord->data[i-1] > dchord->data[(i - 1)-1]) && (dchord->data[i-1] > dchord->data[(i + 1)-1]) && (dchord->data[i-1] > max(rd - 0.0005, reps)) )
        {
            superknots2->data[i - 1] = t->data[i - 1];
        }
    }

    //superknots2(superknots2 == 0) = [];
    m = 0;
    superknots1Data = (double*)malloc(superknots2->column * sizeof(double));
    for (i = 0; i < superknots2->column; i++)
    {
        if (fabs(superknots2->data[i]) > 0.00000001)
        {
            superknots1Data[m] = superknots2->data[i];
            m = m + 1;
        }
    }
    n = 0;
    for (i = 0; i < superknots2->column; i++)
    {
        if (fabs(superknots2->data[i]) > 0.00000001)
        {
            superknots2->data[n] = superknots1Data[n];
            n = n + 1;
        }
    }
    superknots2->column = m;
    free(superknots1Data);

    //superknots = unique(sort([superknots1 superknots2]));
    superknots1Tmp  = Matrix_gen1(1, superknots1->column + superknots2->column, 0.00);
    for (i=0; i< superknots1->column + superknots2->column; i++)
    {
        if ( i < superknots1->column )
        {
            superknots1Tmp->data[i] = superknots1->data[i];
        }
        else
        {
            superknots1Tmp->data[i] = superknots2->data[i - superknots1->column];
        }
    }

    ZXSort(superknots1Tmp->data, superknots1Tmp->column + 1 );
    i = ZXUnique(superknots1Tmp->data, superknots1Tmp->column );
    superknots1Tmp->column = i;
    
    //superknots = superknots1Tmp;
    //Primsuperknots = Matrix_gen1( 1, superknots1Tmp->column - 1, 0.00);
    //% 防止新添加的节点互相离得太近
    /*
    while length(find(diff(superknots) < 0.003)) >= 1
            kk = find(diff(superknots) < 0.003, 1, 'first');
            if kk == length(superknots) - 1
                superknots = superknots(1:kk);
            else
                superknots = [superknots(1:kk) superknots(kk + 2:end)];
            end
    end
    */
    Primsuperknots->column = superknots1Tmp->column - 1;
    ZXdiff(superknots1Tmp->data, Primsuperknots->data, superknots1Tmp->column );
    i = 0;
    for (m = 0; m < Primsuperknots->column; m++)
    {
        if (1 == Primsuperknots->row)
        {
             if (Primsuperknots->data[m] < 0.003)
             {
                i = i + 1;
             }
        }
        /*
        else
        {
            for (n = 0; n < superknots->row; n++)
            {
               if (superknots->data[m * superknots->row + n] < 0.003)
               {
                    i = m * superknots->row + n + 1;
               }
            }
        }
        */
    }

    /*
    if kk == length(superknots) - 1
        superknots = superknots(1:kk);
    else
        superknots = [superknots(1:kk) superknots(kk + 2:end)];
    end
    */
    
    while ( i >= 1 ) 
    {
        m  = 0;
        kk = 0;
        for (m=0;  m < Primsuperknots->column; m++)
        {
            if (Primsuperknots->data[m] < 0.003 )
            {
                kk = m + 1;
                break;
            }
            m = m + 1;
        }

        if ( kk == superknots1Tmp->column - 1)
        {
            superknots1Tmp->data[0] = superknots1Tmp->data[kk-1];
        }
        else
        {
            superknots1Tmp->data[0] = superknots1Tmp->data[kk-1];
            superknots1Tmp->data[1] = superknots1Tmp->data[((superknots1Tmp->column-1) * superknots1Tmp->row) + kk+2 - 1 ];
            superknots1Tmp->column  = 2;
        }
    }

     //% 防止新添加的节点离原节点太近
    //kk = zeros(1, length(superknots));
    Matkk = Matrix_gen1(1, superknots1Tmp->column, 0.00);
    /*
    for i = 1:length(superknots)
        if sum(abs(superknots(i) - sp.knots) < 0.003) > 0
            kk(i) = i;
        end
    end
    */
    for (i=0; i< superknots1Tmp->column; i++)
    {
        double AbsSum = 0.00;
        m = 0;
        for ( j=0; j <sp->knots->column; j++)
        {
            if (fabs(superknots1Tmp->data[i] - sp->knots->data[j]) < 0.003)
            {
                m = m + 1;
            }
        }
        if ( m > 0 )
        {
            Matkk->data[i] = i+1;
        }
    }
    //kk(kk == 0) = [];
    iMatkkZeroCnt = 0;
    for ( i=0; i < Matkk->column; i++ )
    {
        if ( fabs(Matkk->data[i]) > 0.00000001 )
        {
            iMatkkZeroCnt = iMatkkZeroCnt + 1;
        }
    }
    Matrix  *MatkkBak;
    MatkkBak      = Matrix_gen1(1, iMatkkZeroCnt, 0.00);
    iMatkkZeroCnt = 0;
    for ( i=0; i < Matkk->column; i++ )
    {
        if ( fabs(Matkk->data[i]) > 0.00000001 )
        {
            MatkkBak->data[iMatkkZeroCnt] = Matkk->data[i];
            iMatkkZeroCnt = iMatkkZeroCnt + 1;
        }
    }

    //superknots = setdiff(superknots, superknots(kk))
    m = 0;
    int IsExist;
    Matrix* SetDiffExist;
    SetDiffExist = Matrix_gen1(1, Matkk->column, 0.00);

    IsExist = 0;
    m       = 0;
    for ( i=0; i< superknots1Tmp->column; i++ )
    {
        if ( (superknots1Tmp->data[i]) != (superknots1Tmp->data[(int)(MatkkBak->data[0]-1)])  )
        {
            SetDiffExist->data[i- IsExist] = superknots1Tmp->data[i];
            m = m + 1;
        }
        else
        {
            IsExist = IsExist + 1;
        }
    }
    SetDiffExist->column = m;

    /*
    if isempty(superknots)
        nsp = sp;
    else
        nknot = sort([sp.knots superknots]);
        ncontr = Binsert(sp, superknots);
        nsp = spmak2(nknot, ncontr);
    end
    */

    if (SetDiffExist->column > 0)
    {
        //nknot = sort([sp.knots superknots]);
        Matrix* Retnknot;
        Retnknot = Matrix_gen1(1, sp->knots->column + SetDiffExist->column, 0.00);
        for (i = 0; i < sp->knots->column + SetDiffExist->column; i++)
        {
            if (i < sp->knots->column)
            {
                Retnknot->data[i] = sp->knots->data[i];
            }
            else
            {
                Retnknot->data[i] = SetDiffExist->data[i - sp->knots->column];
            }
        }

        ZXSort(Retnknot->data, Retnknot->column + 1);
        //ncontr = Binsert(sp, superknots);
        ncontr = Matrix_gen1(5, sp->number + SetDiffExist->column, 0.00);
        Binsert(ncontr, sp, SetDiffExist);
        //nsp = spmak2(nknot, ncontr);
        spmak2(Retnknot, ncontr, nsp);

    }
    else
    {
        nsp = sp;
    }
    Primsuperknots->column = SetDiffExist->column;
    for (i = 0; i < Primsuperknots->column; i++)
    {
        Primsuperknots->data[i] = SetDiffExist->data[i];
    }
    M_free(superknots1);
    M_free(superknots2);
    M_free(superknots1Tmp);
    M_free(Matkk);
    M_free(MatkkBak);
    M_free(SetDiffExist);
}
