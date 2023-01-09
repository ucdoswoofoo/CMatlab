#include "bangc.h"
#include "ZXMat.h"
#include "matrix.h"

/*
% sp原样条曲线，nt由运动约束放缩节点向量后的新节点向量，alpha, beta每个节点区间的放缩系数
% 输入：1 sp 样条结构体 2 nt double型数组新节点 3 alpha double型数组  4 beta double型数组
% 输出：ncoef使新曲线与原曲线形状相同的对应的控制点
function ncoef = bangc(sp, nt, alpha, beta)
*/
void  bangc(Matrix *ncoef, const  stusp1 sp, const Matrix nt, const Matrix alpha, const Matrix beta)
{
	int n;
	Matrix *t, *coef;

	n    = 0;
	n    = sp.number;
	t    = sp.knots;
	coef = sp.coefs;

	//ncoef = zeros(5, n);
	//ncoef = Matrix_gen1(5, n, 0.00);
	//for i = 4:n
	int  i, j;
	i = 0;
	j = 0;
	for (i = 4; i <= n; i++)
	{
		double g1, g2, g3, g4, ng1, ng2, ng3, ng4;

		//g1 = (t(i + 3) - t(i)) * (t(i + 2) - t(i)) * (t(i + 1) - t(i));
		g1   = (t->data[i + 3 - 1] - t->data[i - 1]) * (t->data[i + 2 - 1] - t->data[i - 1]) * (t->data[i+1 - 1] - t->data[i - 1]);
		//g2 = (t(i + 2) - t(i - 1)) * (t(i + 1) - t(i - 1)) * (t(i + 1) - t(i));
		g2   = (t->data[i+2 - 1]-t->data[i-1 - 1]) * ( t->data[i+1 - 1] -t->data[i-1 - 1]) * ( t->data[i+1 - 1] - t->data[i - 1]);
		//g3 = (t(i + 2) - t(i - 1)) * (t(i + 2) - t(i)) * (t(i + 1) - t(i));
		g3   = (t->data[i+2 - 1] - t->data[i-1 - 1]) * (t->data[i+2 - 1]-t->data[i - 1]) * ( t->data[i+1 - 1] - t->data[i - 1]);
		//g4 = (t(i + 1) - t(i - 2)) * (t(i + 1) - t(i - 1)) * (t(i + 1) - t(i));
		g4   = (t->data[i+1 - 1] - t->data[i-2 - 1]) * (t->data[i+1 - 1] - t->data[i-1 - 1]) * (t->data[i+1 - 1] - t->data[i - 1]);
		

		//ng1 = (nt(i + 3) - nt(i)) * (nt(i + 2) - nt(i)) * (nt(i + 1) - nt(i));
		ng1   = (nt.data[i+3 - 1] - nt.data[i - 1]) * ( nt.data[i+2 - 1]-nt.data[i - 1]) * ( nt.data[i+1 - 1] - nt.data[i - 1]);
		//ng2 = (nt(i + 2) - nt(i - 1)) * (nt(i + 1) - nt(i - 1)) * (nt(i + 1) - nt(i));
		ng2   = (nt.data[i+2-1] - nt.data[i-1-1]) * (nt.data[i+1-1] - nt.data[i-1-1]) * (nt.data[i+1-1] - nt.data[i-1]);
		//ng3 = (nt(i + 2) - nt(i - 1)) * (nt(i + 2) - nt(i)) * (nt(i + 1) - nt(i));
		ng3   = (nt.data[i+2-1] - nt.data[i-1-1]) * (nt.data[i+2-1] - nt.data[i-1]) * (nt.data[i+1-1] - nt.data[i-1]);
		//ng4 = (nt(i + 1) - nt(i - 2)) * (nt(i + 1) - nt(i - 1)) * (nt(i + 1) - nt(i));
		ng4   = (nt.data[i+1-1] - nt.data[i-2-1]) * (nt.data[i+1-1] - nt.data[i-1-1]) * (nt.data[i+1-1] - nt.data[i-1]);

		//% 新控制点 = A * B * (C的逆)
		Matrix* A;		
		A = M_Cut(coef, 1, coef->row, i-3, i);

		double c24, c34;

		c24 = (-1) * pow(nt.data[i - 1 - 1], 3) / ng2 - pow(nt.data[i + 2 - 1], 3) / ng3 - pow(nt.data[i + 1 - 1], 3) / ng4 + (nt.data[i - 1 - 1] + nt.data[i - 1] + nt.data[i + 2 - 1]) / (nt.data[i + 1 - 1] - nt.data[i - 1]) + 1;
		c34 = pow(nt.data[i - 1], 3) / ng1 + pow(nt.data[i - 1 - 1], 3) / ng2 + pow(nt.data[i + 2 - 1], 3) / ng3 - ( nt.data[i-1-1] + nt.data[i-1] + nt.data[i+2-1]) / (nt.data[i + 1 - 1] - nt.data[i - 1]);

		Matrix* C;
		C = Matrix_gen1(4, 4, 0.00);
		/*
		C = [-1 / ng4 3 * nt(i + 1) / ng4 - 3 * nt(i + 1) ^ 2 / ng4 nt(i + 1) ^ 3 / ng4;
		1 / ng2 + 1 / ng3 + 1 / ng4 - 3 * (nt(i - 1) / ng2 + nt(i + 2) / ng3 + nt(i + 1) / ng4)  3 * (nt(i - 1) ^ 2 / ng2 + nt(i + 2) ^ 2 / ng3 + nt(i + 1) ^ 2 / ng4 - 1 / (nt(i + 1) - nt(i))) c24;
		-1 / ng1 - 1 / ng2 - 1 / ng3 3 * (nt(i) / ng1 + nt(i - 1) / ng2 + nt(i + 2) / ng3) 3 * (-nt(i) ^ 2 / ng1 - nt(i - 1) ^ 2 / ng2 - nt(i + 2) ^ 2 / ng3 + 1 / (nt(i + 1) - nt(i))) c34;
		1 / ng1 - 3 * nt(i) / ng1 3 * nt(i) ^ 2 / ng1 - nt(i) ^ 3 / ng1];
		*/

		C->data[0]  = -1/ng4;
		C->data[4]  = 3 * nt.data[i+1-1]/ng4;
		C->data[8]  = -3*pow(nt.data[i+1-1], 2)/ng4 ;
		C->data[12] = pow(nt.data[i+1-1], 3) / ng4;
		C->data[1]  = 1 / ng2 + 1 / ng3 + 1 / ng4;
		C->data[5]  = -3 * (nt.data[i -1- 1] / ng2 + nt.data[i + 2-1] / ng3 + nt.data[i + 1-1] / ng4);
		C->data[9]  = 3 * (pow(nt.data[i - 1-1], 2) / ng2 + pow(nt.data[i + 2-1], 2) / ng3 + pow(nt.data[i + 1-1], 2) / ng4 - 1 / (nt.data[i + 1-1] - nt.data[i-1]));
		C->data[13] = c24;
		C->data[2]  = -1 / ng1 - 1 / ng2 - 1 / ng3;
		C->data[6]  = 3 * (nt.data[i-1] / ng1 + nt.data[i - 1-1] / ng2 + nt.data[i + 2-1] / ng3);
		C->data[10] = 3 * (-pow(nt.data[i-1] , 2) / ng1 - pow(nt.data[i - 1-1],  2) / ng2 - pow(nt.data[i + 2-1] , 2) / ng3 + 1 / (nt.data[i + 1 - 1] - nt.data[i -1]));
		C->data[14] = c34;
		C->data[3]  = 1 / ng1;
		C->data[7]  = -3 * nt.data[i-1] / ng1;
		C->data[11] = 3 * pow(nt.data[i-1], 2) / ng1;
		C->data[15] = -pow(nt.data[i-1], 3) / ng1;

		double b24, b34;

		b24 = -pow(t->data[i -1- 1], 3) / g2 - pow(t->data[i + 2-1], 3) / g3 - pow(t->data[i + 1-1], 3) / g4 + (t->data[i -1- 1] + t->data[i-1] + t->data[i + 2-1]) / (t->data[i + 1-1] - t->data[i-1]) + 1;
		b34 = pow(t->data[i-1], 3) / g1 + pow(t->data[i -1- 1], 3) / g2 + pow(t->data[i + 2-1], 3) / g3 - (t->data[i -1- 1] + t->data[i-1] + t->data[i + 2-1]) / (t->data[i + 1-1] - t->data[i-1]);
		/*
		b24 = -t.CellData[(i - 1) % (t.iRow)][(i - 1) / (t.iRow)] * t.CellData[(i - 1) % (t.iRow)][(i - 1) / (t.iRow)] * t.CellData[(i - 1) % (t.iRow)][(i - 1) / (t.iRow)] / g2 - t.CellData[(i + 2) % (t.iRow)][(i + 2) / (t.iRow)] * t.CellData[(i + 2) % (t.iRow)][(i + 2) / (t.iRow)] * t.CellData[(i + 2) % (t.iRow)][(i + 2) / (t.iRow)] / g3 - t.CellData[(i + 1) % (t.iRow)][(i + 1) / (t.iRow)] * t.CellData[(i + 1) % (t.iRow)][(i + 1) / (t.iRow)] * t.CellData[(i + 1) % (t.iRow)][(i + 1) / (t.iRow)] / g4 + (t.CellData[(i - 1) % (t.iRow)][(i - 1) / (t.iRow)] + t.CellData[(i) % (t.iRow)][(i) / (t.iRow)] + t.CellData[(i + 2) % (t.iRow)][(i + 2) / (t.iRow)]) / (t.CellData[(i + 1) % (t.iRow)][(i + 1) / (t.iRow)] - t.CellData[(i) % (t.iRow)][(i) / (t.iRow)]) + 1;
		b34 = t.CellData[(i) % (t.iRow)][(i) / (t.iRow)] * t.CellData[(i) % (t.iRow)][(i) / (t.iRow)] * t.CellData[(i) % (t.iRow)][(i) / (t.iRow)] / g1 + t.CellData[(i - 1) % (t.iRow)][(i - 1) / (t.iRow)] * t.CellData[(i - 1) % (t.iRow)][(i - 1) / (t.iRow)] * t.CellData[(i - 1) % (t.iRow)][(i - 1) / (t.iRow)] / g2 + t.CellData[(i + 2) % (t.iRow)][(i + 2) / (t.iRow)] * t.CellData[(i + 2) % (t.iRow)][(i + 2) / (t.iRow)] * t.CellData[(i + 2) % (t.iRow)][(i + 2) / (t.iRow)] / g3 - (t.CellData[(i - 1) % (t.iRow)][(i - 1) / (t.iRow)] + t.CellData[(i) % (t.iRow)][(i) / (t.iRow)] + t.CellData[(i + 2) % (t.iRow)][(i + 2) / (t.iRow)]) / (t.CellData[(i + 1) % (t.iRow)][(i + 1) / (t.iRow)] - t.CellData[(i) % (t.iRow)][(i) / (t.iRow)]);
		*/
		Matrix* Btemp;

		Btemp           = Matrix_gen1(4, 4, 0.00);
		Btemp->data[0]  = -1 / g4;
		Btemp->data[4]  = 3 * t->data[i + 1-1] / g4;
		Btemp->data[8]  = -3 * pow(t->data[i + 1-1], 2) / g4;
		Btemp->data[12] = pow(t->data[i + 1-1], 3) / g4;
		Btemp->data[1]  = 1 / g2 + 1 / g3 + 1 / g4;
		Btemp->data[5]  = -3 * (t->data[i - 1-1] / g2 + t->data[i + 2-1] / g3 + t->data[i + 1-1] / g4);
		Btemp->data[9]  = 3 *( pow(t->data[i - 1-1], 2) / g2 + pow(t->data[i + 2-1], 2) / g3 + pow(t->data[i + 1-1], 2) / g4 - (1 / (t->data[i + 1-1] - t->data[i-1])));
		Btemp->data[13] = b24;
		Btemp->data[2]  = -1 / g1 - 1 / g2 - 1 / g3;
		Btemp->data[6]  = 3 * (t->data[i-1] / g1 + t->data[i - 1-1] / g2 + t->data[i + 2-1] / g3);
		Btemp->data[10] = 3 * (-pow(t->data[i-1], 2) / g1 - pow(t->data[i - 1-1], 2) / g2 - pow(t->data[i + 2-1], 2) / g3 + 1 / (t->data[i + 1-1] - t->data[i-1]));
		Btemp->data[14] = b34;
		Btemp->data[3]  = 1 / g1;
		Btemp->data[7]  = -3 * t->data[i-1] / g1;
		Btemp->data[11] = 3 * pow(t->data[i-1], 2) / g1;
		Btemp->data[15] = -pow(t->data[i-1], 3) / g1;

		/*
		Btemp = [-1 / g4 3 * t(i + 1) / g4 - 3 * t(i + 1) ^ 2 / g4 t(i + 1) ^ 3 / g4;
		1 / g2 + 1 / g3 + 1 / g4 - 3 * (t(i - 1) / g2 + t(i + 2) / g3 + t(i + 1) / g4)  3 * (t(i - 1) ^ 2 / g2 + t(i + 2) ^ 2 / g3 + t(i + 1) ^ 2 / g4 - 1 / (t(i + 1) - t(i))) b24;
		-1 / g1 - 1 / g2 - 1 / g3 3 * (t(i) / g1 + t(i - 1) / g2 + t(i + 2) / g3) 3 * (-t(i) ^ 2 / g1 - t(i - 1) ^ 2 / g2 - t(i + 2) ^ 2 / g3 + 1 / (t(i + 1) - t(i))) b34;
		1 / g1 - 3 * t(i) / g1 3 * t(i) ^ 2 / g1 - t(i) ^ 3 / g1];
		*/

		Matrix* B;
		
		B = Matrix_gen1(4, 4, 0.00);
		B->data[0] = Btemp->data[0] * pow(alpha.data[i - 3-1], 3);
		B->data[1] = Btemp->data[1] * pow(alpha.data[i - 3 - 1], 3);
		B->data[2] = Btemp->data[2] * pow(alpha.data[i - 3 - 1], 3);
		B->data[3] = Btemp->data[3] * pow(alpha.data[i - 3 - 1], 3);
		B->data[4] = Btemp->data[0] * 3 * pow(alpha.data[i - 3 - 1], 2) * beta.data[i - 3 - 1] + Btemp->data[4] * pow(alpha.data[i - 3 - 1], 2);
		B->data[5] = Btemp->data[1] * 3 * pow(alpha.data[i - 3 - 1], 2) * beta.data[i - 3 - 1] + Btemp->data[5] * pow(alpha.data[i - 3 - 1], 2);
		B->data[6] = Btemp->data[2] * 3 * pow(alpha.data[i - 3 - 1], 2) * beta.data[i - 3 - 1] + Btemp->data[6] * pow(alpha.data[i - 3 - 1], 2);
		B->data[7] = Btemp->data[3] * 3 * pow(alpha.data[i - 3 - 1], 2) * beta.data[i - 3 - 1] + Btemp->data[7] * pow(alpha.data[i - 3 - 1], 2);
		B->data[8] = Btemp->data[0] * 3 * alpha.data[i - 3-1] * pow(beta.data[i - 3-1], 2) + Btemp->data[4] * 2 * alpha.data[i - 3-1] * beta.data[i - 3 - 1] + Btemp->data[8] * alpha.data[i - 3 - 1];
		B->data[9] = Btemp->data[1] * 3 * alpha.data[i - 3 - 1] * pow(beta.data[i - 3 - 1], 2) + Btemp->data[5] * 2 * alpha.data[i - 3 - 1] * beta.data[i - 3 - 1] + Btemp->data[9] * alpha.data[i - 3 - 1];
		B->data[10] = Btemp->data[2] * 3 * alpha.data[i - 3 - 1] * pow(beta.data[i - 3 - 1], 2) + Btemp->data[6] * 2 * alpha.data[i - 3 - 1] * beta.data[i - 3 - 1] + Btemp->data[10] * alpha.data[i - 3 - 1];
		B->data[11] = Btemp->data[3] * 3 * alpha.data[i - 3 - 1] * pow(beta.data[i - 3 - 1], 2) + Btemp->data[7] * 2 * alpha.data[i - 3 - 1] * beta.data[i - 3 - 1] + Btemp->data[11] * alpha.data[i - 3 - 1];
		B->data[12] = Btemp->data[0] * pow(beta.data[i - 3 - 1], 3) + Btemp->data[4] * pow(beta.data[i - 3-1], 2) + Btemp->data[8] * beta.data[i - 3 - 1] + Btemp->data[12];
		B->data[13] = Btemp->data[1] * pow(beta.data[i - 3 - 1], 3) + Btemp->data[5] * pow(beta.data[i - 3 - 1], 2) + Btemp->data[9] * beta.data[i - 3 - 1] + Btemp->data[13];
		B->data[14] = Btemp->data[2] * pow(beta.data[i - 3 - 1], 3) + Btemp->data[6] * pow(beta.data[i - 3 - 1], 2) + Btemp->data[10] * beta.data[i - 3 - 1] + Btemp->data[14];
		B->data[15] = Btemp->data[3] * pow(beta.data[i - 3-1], 3) + Btemp->data[7] * pow(beta.data[i - 3 - 1], 2) + Btemp->data[11] * beta.data[i - 3-1] + Btemp->data[15];

		//D = A * B * C ^ (-1);
		Matrix* D = M_Inverse(C);
		Matrix* Mul1 = M_mul(A, B);
		Matrix* Mul2 = M_mul(Mul1, D);

		//ncoef(:, i - 3 : i) = D;
		for (int n = (i - 3); n <= i; n++)
		{
			for (int m = 0; m < ncoef->row; m++)
			{
				ncoef->data[(n - (i-3)) * Mul2->row + m] = Mul2->data[(n-(i-3)) * Mul2->row + m];
			}
		}
		
		M_free(A);
		M_free(C);
		M_free(D);
		M_free(Mul1);
		M_free(Mul2);
		M_free(B);
		M_free(Btemp);
	}

}


