#include "spremove2.h"
#include "ZXMat.h"
#include "matrix.h"
#include "spmak2.h"

void  spremove2(stusp1* spn, const stusp1* sp, const int r)
{
	Matrix * knot, * coef, * temp1, * temp2, * knotn, * coefn;
	int p, i, j, m;
	double u, alf1, alf2;

	alf1 = 0.00;
	alf2 = 0.00;
	i = 0;
	j = 0;
	m = 0;
	knot = sp->knots;
	coef = sp->coefs;

	//p = 3;
	//u = knot(r);
	p = 3;
	u = knot->data[r-1];

	//% alphar - pºÍalphar - p + 1
	//alf1 = (u - knot(r - p)) / (knot(r + 1) - knot(r - p));
	//alf2 = (u - knot(r - p + 1)) / (knot(r + 2) - knot(r - p + 1));
	alf1 = (u - knot->data[r - p - 1]) / (knot->data[(r + 1)-1] - knot->data[(r - p)-1]);
	alf2 = (u - knot->data[(r - p + 1)-1]) / (knot->data[(r + 2)-1] - knot->data[(r - p + 1)-1]);


	//temp1 = (coef(:, r - p) - (1 - alf1) * coef(:, r - p - 1)) / alf1;
	//temp2 = (coef(:, r - p + 1) - (1 - alf2) * temp1) / alf2;
	temp1 = Matrix_gen1(coef->row, 1, 0.00);
	temp2 = Matrix_gen1(coef->row, 1, 0.00);
	for (j = 0; j < coef->row; j++)
	{
		temp1->data[j] = (coef->data[coef->row * (r - p-1) + j] - (1 - alf1) * coef->data[coef->row * (r - p - 1-1) + j]) / alf1;
	}
	for (j = 0; j < coef->row; j++)
	{
		temp2->data[j] = (coef->data[coef->row * (r - p + 1-1) + j] - (1 - alf2) * temp1->data[j]) / alf2;
	}

	//knotn = [knot(1:r - 1) knot(r + 1:end)];
	//m = 0;
	knotn = Matrix_gen1(knot->row, knot->column-1, 0.00);
	for (i = 0; i < r-1; i++)
	{
		knotn->data[i] = knot->data[i];
		//m = i + 1;
	}
	for (i=r+1-1; i < knot->column;  i++)
	{
		knotn->data[i - 1] = knot->data[i];
	}
	//coefn = [coef(:, 1 : r - p - 1) temp1 temp2 coef(:, r - p + 3 : end)];
	coefn = Matrix_gen1(coef->row, (r-p-1)-(1-1)+temp1->column+temp2->column + coef->column - (r-p+3-1), 0.00);
	m = 0;
	for (i = 0; i < r-p-1; i++)
	{
		for (j=0; j < coef->row; j++)
		{
			coefn->data[i*coef->row + j] = coef->data[i * coef->row + j];
		}
	}
	m = (r-p-1) * coefn->row;
	for (i = 0; i < temp1->column; i++)
	{
		for (j = 0; j < temp1->row; j++)
		{
			coefn->data[m + i * temp1->row + j] = temp1->data[i * temp1->row + j];
		}
	}
	m = m + temp1->column * temp1->row;
	for (i = 0; i < temp2->column; i++)
	{
		for (j = 0; j < temp2->row; j++)
		{
			coefn->data[m + i * temp2->row + j] = temp2->data[i * temp2->row + j];
		}
	}
	m = m + temp2->column * temp2->row;
	for (i = r - p + 3-1; i < coef->column; i++)
	{
		for (j = 0; j < coef->row; j++)
		{
			coefn->data[ m + (i - (r - p + 3-1)) * coef->row + j] = coef->data[i * coef->row + j];
		}
	}

	M_free(temp1);
	M_free(temp2);
	//spn = spmak2(knotn, coefn);
	spmak2(knotn, coefn, spn);
	//spn->number = coefn->column;
	//spn->knots = knotn;
	//spn->coefs = coefn;

}


