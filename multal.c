#include "multal.h"

/*
% 放大alpha，相当于新的alpha向量是所有元素都是1.几
% 输入：1 alpha double型1 * 1    2 sp 样条结构体
function[tknot, nalpha, nbeta] = multal(alpha, sp)
*/
void  multal(Matrix* tknot, Matrix* nalpha, Matrix* nbeta, const double alpha, const stusp1* sp)
{
	//n = sp.number;
	int n, i;
	Matrix* knot;

	i    = 0;
	n    = sp->number;
	knot = sp->knots;

	/*
	for i = 5:n
	tknot(i) = (knot(i) - nbeta(i - 4)) / alpha;
	nbeta(i - 3) = knot(i) - alpha * tknot(i);
	end
	*/
	
	for (i=5; i<=n; i++)
	{
		if ((tknot->column >= i - 1) && (knot->column >= i - 1) && (nbeta->column >= i - 5))
		{
			tknot->data[i - 1] = (knot->data[i - 1] - nbeta->data[i - 5]) / alpha;
		}
		
		if ( (nbeta->column >= i-4) & (knot->column >= i-1) && (tknot->column >= i-1) )
		{
			nbeta->data[i - 4] = knot->data[i - 1] - alpha * tknot->data[i - 1];
		}	
	}
	
	/*
	tknot(n + 1) = (knot(n + 1) - nbeta(n - 3)) / alpha;
	tknot(n + 2:end) = tknot(n + 1);
	*/
	if ( (tknot->column >= n) && (knot->column >= n) && (nbeta->column >= n-4) )
	{
		tknot->data[n] = (knot->data[n] - nbeta->data[n - 4]) / alpha;
	}
	
	for (i=n+2; i<= tknot->column; i++)
	{
		if ( ( tknot->column > i-1) && (tknot->column > n) )
		tknot->data[i - 1] = tknot->data[n];
	}
}


