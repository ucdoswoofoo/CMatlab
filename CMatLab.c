// CMatLab.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include "stdlib.h"
#include "stdio.h"
#include "ZXMat.h"
#include "bangc.h"
#include "matrix.h"
#include "spmak2.h"
#include "spremove2.h"
#include "globaltrans.h"
#include "addknot.h"
#include "Binsert.h"
#include "sys\types.h"
#include <time.h>

//#define ALPHA 3.95


int main()
{


	float n;
	int ch, rani;

	//printf("\nRandom numbers between 0 to 1\n");
	/*
	for (rani =0; rani <25; rani++)
	{
		n = (float)rand() / RAND_MAX;
		printf("Random number[%d]: %f\n", rani, n);
		//printf("\nDo you want to generate again (1/0): ");
		//scanf("%d", &ch);
	} //while (ch == 1);

	system("pause");
	exit(0);
	*/


	//开始测试 addknot
	int i, iCol, iRow;
	i    = 0;
	iCol = 0;
	iRow = 0;
	double pii  = 3.141592653589793;
	double reps = 0.03 * pii / 180;
	double eps  = 0.005;

	Matrix* knot, * coef, * pkd, * pkpar, * dchord;
	knot = Matrix_gen1(1, 9, 0.00);
	knot->data[0] = 0;
	knot->data[1] = 0;
	knot->data[2] = 0;
	knot->data[3] = 0;
	knot->data[4] = 0.5;
	knot->data[5] = 1;
	knot->data[6] = 1;
	knot->data[7] = 1;
	knot->data[8] = 1;


	/*
	coef = Matrix_gen1(5, 5, 0.00);

	for (i=0; i<25; i++)
	{
		n = (float)rand() / RAND_MAX;
		//printf("Random number: %f\n", n);
		coef->data[i] = n;
	}
	*/
	MATRIX_TYPE _MatTmp1[25] = { 0.162182308193243,0.794284540683907,0.311215042044805,0.528533135506213,0.165648729499781,
								 0.601981941401637,0.262971284540144,0.654079098476782,0.689214503140008,0.748151592823710,
								 0.450541598502498,0.0838213779969326,0.228976968716819,0.913337361501670,0.152378018969223,
								 0.825816977489547,0.538342435260057,0.996134716626886,0.0781755287531837,0.442678269775446,
								 0.106652770180584,0.961898080855054,0.00463422413406744,0.774910464711502,0.817303220653433  };

	coef = Matrix_gen(5, 5, _MatTmp1);
	for (rani = 0; rani < 25; rani++)
	{
		n = (float)rand() / RAND_MAX;
		coef->data[rani] = n;
		//printf("Random number[%d]: %f\n", rani, n);
		//printf("\nDo you want to generate again (1/0): ");
		//scanf("%d", &ch);
	}
	stusp1* spaddknot;
	spaddknot = (stusp1*)malloc(sizeof(stusp1));
	spaddknot->number = 5;
	spaddknot->coefs = coef;
	spaddknot->knots = knot;
	pkd = Matrix_gen1(1, 4, 0.00);
	pkd->data[0] = 0.01;
	pkd->data[1] = 0.02;
	pkd->data[2] = 0.036;
	pkd->data[3] = 0.04;
	double wchorderr = 0.04;
	double wdataerr = 0.01;
	pkpar = Matrix_gen1(1, 4, 0.00);
	pkpar->data[0] = 0.3;
	pkpar->data[1] = 0.4;
	pkpar->data[2] = 0.5003;
	pkpar->data[3] = 0.7;
	dchord = Matrix_gen1(1, 5, 0.00);
	dchord->data[0] = 0.002;
	dchord->data[1] = 0.003;
	dchord->data[2] = 0.001;
	dchord->data[3] = 0.004;
	dchord->data[4] = 0.0003;
	double rd;
	rd = 0.004;

	stusp1* addknotnsp;
	Matrix* superknots;

	addknotnsp = (stusp1*)malloc(sizeof(stusp1));


	/*
	
	knot = [0 0 0 0 0.5 1 1 1 1];
	coef = rand(5, 5);
	sp = spmak2(knot, coef);
	pkd = [0.01 0.02 0.036 0.04];
	double wchorderr = 0.04;
	double wdataerr = 0.01;
	pkpar = [0.3 0.4 0.5003 0.7];
	dchord = [0.002 0.003 0.001 0.004 0.0003];
	rd = 0.004;
	[nsp, superknots] = addknot(sp, pkd, wchorderr, wdataerr, eps, pkpar, dchord, rd, reps);
	*/
	printf("Begin Call addknot .................................................\n");
	printf("wchorderr=%lf, wdataerr=%lf, eps=%lf, rd=%lf, reps=%lf\n", wchorderr, wdataerr, eps, rd, reps);
	printf("spaddknot->number=[%d]\n", spaddknot->number);
	for (iCol = 0; iCol < spaddknot->knots->column; iCol++)
	{
		printf("spaddknot->knots->data[%d] = %.4lf\n", iCol, spaddknot->knots->data[iCol]);
	}
	for (iCol = 0; iCol < spaddknot->coefs->column; iCol++)
	{
		for (iRow=0; iRow< spaddknot->coefs->row; iRow++)
		{
			printf("spaddknot->coefs->data[Row=%d][Column=%d] = %.19lf\n", iRow, iCol, spaddknot->coefs->data[iCol * spaddknot->coefs->row + iRow]);
		}
	}
	for (iCol = 0; iCol < pkd->column; iCol++)
	{
		printf("pkd->data[%d] = %.4lf\n", iCol, pkd->data[iCol]);
	}
	for (i = 0; i < pkpar->column; i++)
	{
		printf("pkpar->data[%d] = %.4lf\n", i, pkpar->data[i]);
	}
	for (i = 0; i < dchord->column; i++)
	{
		printf("dchord->data[%d] = %.4lf\n", i, dchord->data[i]);
	}

	superknots = Matrix_gen1(1, pkd->column + 1000, 0);

	printf("Call addknot  ***************************************************\n");
	addknot(addknotnsp, superknots, spaddknot, pkd, wchorderr, wdataerr, eps, pkpar, dchord, rd, reps);

	printf("addknot return value : **********************************************\n");
	
	
	for (i = 0; i < superknots->column; i++)
	{
		printf("superknots->data[%d] = %.4lf\n", i, superknots->data[i]);
	}
	
	printf("addknotnsp->number[%d]\n",  addknotnsp->number);
	for (i = 0; i < addknotnsp->coefs->column; i++)
	{
		for (iRow = 0; iRow < addknotnsp->coefs->row; iRow++)
		{
			printf("addknotnsp->coefs->data[Row=%d][Column=%d] = %.19lf\n", iRow, i, addknotnsp->coefs->data[i * addknotnsp->coefs->row + iRow]);
		}
	}
	for (i = 0; i < addknotnsp->knots->column; i++)
	{
		printf("addknotnsp->knots->data[%d] = %.19lf\n", i, addknotnsp->knots->data[i]);
	}
	M_free(pkd);
	M_free(pkpar);
	M_free(dchord);
	M_free(spaddknot->coefs);
	M_free(spaddknot->knots);
	free(spaddknot);

	M_free(addknotnsp->coefs);
	M_free(addknotnsp->knots);
	free(addknotnsp);
	M_free(superknots);
	printf("Call addknot Over **************************************************\n");
	//exit(0);
	//system("pause");
	//exit(0);
	//结束测试 addknot

	//开始测试  spremove2
	Matrix* spremove2knot, * spremove2coef;
	spremove2knot = Matrix_gen1(1, 8, 0.00);
	spremove2knot->data[0] = 0;
	spremove2knot->data[1] = 0;
	spremove2knot->data[2] = 0;
	spremove2knot->data[3] = 0;
	spremove2knot->data[4] = 1;
	spremove2knot->data[5] = 1;
	spremove2knot->data[6] = 1;
	spremove2knot->data[7] = 1;

	/*
	spremove2coef = Matrix_gen1(5, 4, 0.00);

	for (i = 0; i < 20; i++)
	{
		n = (float)rand() / RAND_MAX;
		//printf("Random number: %f\n", n);
		spremove2coef->data[i] = n;
	}
	*/

	MATRIX_TYPE _MatTmp[20] = { 0.351659507062997,0.830828627896291,0.585264091152724,0.549723608291140,0.917193663829810,
								0.285839018820374,0.757200229110721,0.753729094278495,0.380445846975357,0.567821640725221,
								0.0758542895630636,0.0539501186666072,0.530797553008973,0.779167230102011,0.934010684229183,
								0.129906208473730,0.568823660872193,0.469390641058206,0.0119020695012414,0.337122644398882 };

	spremove2coef = Matrix_gen(5, 4, _MatTmp);
	for (rani = 0; rani < 20; rani++)
	{
		n = (float)rand() / RAND_MAX;
		spremove2coef->data[rani] = n;
		//printf("Random number[%d]: %f\n", rani, n);
	}
	//M_print(spremove2coef);
	//exit(0);
	//knot = [0 0 0 0 1 1 1 1];
	//coef = rand(5, 4);
	stusp1* spremove2sp, * spremove2nsp, * spremove2nspend;
	spremove2sp = (stusp1*)malloc(sizeof(stusp1));
	//spremove2sp = spmak2(spremove2knot, spremove2coef);
	spmak2(spremove2knot, spremove2coef, spremove2sp);
	//ncontr = Binsert(spremove2sp, 0.5);

	Matrix* spremove2ncontr, * spremove2superKnots;

	spremove2ncontr = Matrix_gen1(5, spremove2sp->number + 1, 0.00);
	spremove2superKnots = Matrix_gen1(1, 1, 0.00);
	spremove2superKnots->data[0] = 0.5;

	Binsert(spremove2ncontr, spremove2sp, spremove2superKnots);

	//spremove2nknotend = [0 0 0 0 0.5 1 1 1 1];
	Matrix* spremove2nknot;
	spremove2nknot = Matrix_gen1(1, 9, 0.00);
	spremove2nknot->data[0] = 0;
	spremove2nknot->data[1] = 0;
	spremove2nknot->data[2] = 0;
	spremove2nknot->data[3] = 0;
	spremove2nknot->data[4] = 0.5;
	spremove2nknot->data[5] = 1;
	spremove2nknot->data[6] = 1;
	spremove2nknot->data[7] = 1;
	spremove2nknot->data[8] = 1;

	//spremove2nsp = spmak2(nknot, ncontr);
	spremove2nsp = (stusp1*)malloc(sizeof(stusp1));
	spmak2(spremove2nknot, spremove2ncontr, spremove2nsp);


	//spremove2nspend = spremove2(spremove2nsp, 5);
	spremove2nspend = (stusp1*)malloc(sizeof(stusp1));
	printf("Begin Call spremove2  ***********************************************************\n");
	printf("spremove2nsp->number=[%d]\n", spremove2nsp->number);
	for (iCol = 0; iCol < spremove2nsp->coefs->column; iCol++)
	{
		for (iRow = 0; iRow < spremove2nsp->coefs->row; iRow++)
		{
			printf("spremove2nsp->coefs->data[Row=%d][Column=%d] = %.19lf\n", iRow, iCol, spremove2nsp->coefs->data[iCol * spremove2nsp->coefs->row + iRow]);
		}
	}
	for (iCol = 0; iCol < spremove2nsp->knots->column; iCol++)
	{
		printf("spremove2nsp->knots->data[%d] = %.19lf\n", iCol, spremove2nsp->knots->data[iCol] );
	}
	printf("Call spremove2  ***************************************************\n");
	spremove2(spremove2nspend, spremove2nsp, 5);
	printf("return spremove2 Over **************************************************************\n");
	printf("spremove2nspend->number=[%d]\n", spremove2nspend->number);
	for (iCol = 0; iCol < spremove2nspend->coefs->column; iCol++)
	{
		for (iRow = 0; iRow < spremove2nspend->coefs->row; iRow++)
		{
			printf("spremove2nspend->coefs->data[Row=%d][Column=%d] = %.19lf\n", iRow, iCol, spremove2nspend->coefs->data[iCol * spremove2nspend->coefs->row + iRow]);
		}
	}
	for (iCol = 0; iCol < spremove2nspend->knots->column; iCol++)
	{
		printf("spremove2nspend->knots->data[Column=%d] = %.19lf\n", iCol, spremove2nspend->knots->data[iCol]);
	}
	printf("Call spremove2 Over .................................................\n");
	printf(".................. .................................................\n");
	//norm(sp.coefs - nsp.coefs)
	Matrix  *coefSub;
	coefSub = M_add_sub( 1, spremove2sp->coefs, 1, spremove2nspend->coefs);

	MATRIX_TYPE   norm;
	
	norm = M_norm(coefSub, 1);

	printf("Matrix  norm = %.19lf\n", norm);

	M_free(spremove2sp->coefs);
	M_free(spremove2sp->knots);
	free(spremove2sp);
	M_free(spremove2nsp->coefs);
	M_free(spremove2nsp->knots);
	free(spremove2nsp);
	M_free(spremove2nspend->coefs);
	M_free(spremove2nspend->knots);
	free(spremove2nspend);
	M_free(coefSub);

	//结束测试   spremove2

	//system("pause");
	//exit(0);

	
	/*
	int ij, ji;
	double a[12] = {0.52, 0.699, 0.53}; // { 11.07431504, 8.824832431, 11.47422964, 9.395644882, 9.395644882, 10.38676979, 10.85168032, 11.28350853, 11.38184929, 11.56150055, 11.64178781, 11.71609794 };
	
	for (ij = 0; ij < 3; ij++)
	{
		printf("%lf\n", a[ij]);
	}
	ZXSort(a, 4);
	printf("sorted\n");

	for (ij = 0; ij < 3; ij++)
	{
		printf("%lf\n", a[ij]);
	}
	printf("Prim Count: %d\n", 12);

	ji = 0;

	ji = ZXUnique( a, 3);
	printf("uniqued:\n");
	for (ij = 0; ij < 12; ij++)
	{
		printf("%lf\n", a[ij]);
	}

	printf("uniqued Count: %d\n", ji);

	system("pause");
	exit(0);
	*/
	
	// std::cout << "Hello World!\n";
	//开始测试  globaltrans
	double pi = 3.14159265;
	double Vmax = 250;
	double Amax = 500;
	double Jmax = 3000;
	double RVmax_A = 300 * pi / 180;
	double RVmax_C = 300 * pi / 180;
	double RAmax = 500 * pi / 180;
	double RJmax = 3000 * pi / 180; 

	//% 误差界
	//double eps = 0.005;

	stusp1* nsp, *sp;

	nsp = (stusp1*)malloc(sizeof(stusp1));
	sp = (stusp1*) malloc(sizeof(stusp1));

	sp->number = 4;
	sp->knots = (Matrix*) malloc(sizeof(Matrix));
	sp->knots->row = 1;
	sp->knots->column = 8;
	sp->knots->data = (double *)malloc( 8 * sizeof(double));
	sp->knots->data[0] = 0;
	sp->knots->data[1] = 0;
	sp->knots->data[2] = 0;
	sp->knots->data[3] = 0;
	sp->knots->data[4] = 1;
	sp->knots->data[5] = 1;
	sp->knots->data[6] = 1;
	sp->knots->data[7] = 1;

	sp->coefs = (Matrix*)malloc(sizeof(Matrix));
	sp->coefs->row = 5;
	sp->coefs->column = 4;
	sp->coefs->data = (double*)malloc( (5 * 4) * sizeof(double));
	/*
	sp->coefs->data[0] = 0.963088539286913;
	sp->coefs->data[1] = 0.546805718738968;
	sp->coefs->data[2] = 0.521135830804002;
	sp->coefs->data[3] = 0.231594386708524;
	sp->coefs->data[4] = 0.488897743920167;
	sp->coefs->data[5] = 0.624060088173690;
	sp->coefs->data[6] = 0.679135540865748;
	sp->coefs->data[7] = 0.395515215668593;
	sp->coefs->data[8] = 0.367436648544477;
	sp->coefs->data[9] = 0.987982003161633;
	sp->coefs->data[10] = 0.0377388662395521;
	sp->coefs->data[11] = 0.885168008202475;
	sp->coefs->data[12] = 0.913286827639239;
	sp->coefs->data[13] = 0.796183873585212;
	sp->coefs->data[14] = 0.0987122786555743;
	sp->coefs->data[15] = 0.261871183870716;
	sp->coefs->data[16] = 0.335356839962797;
	sp->coefs->data[17] = 0.679727951377338;
	sp->coefs->data[18] = 0.136553137355370;
	sp->coefs->data[19] = 0.721227498581740;
	*/
	sp->coefs->data[0] = 0.162182308193243;
	sp->coefs->data[1] = 0.794284540683907;
	sp->coefs->data[2] = 0.311215042044805;
	sp->coefs->data[3] = 0.528533135506213;
	sp->coefs->data[4] = 0.165648729499781;
	sp->coefs->data[5] = 0.601981941401637;
	sp->coefs->data[6] = 0.262971284540144;
	sp->coefs->data[7] = 0.654079098476782;
	sp->coefs->data[8] = 0.689214503140008;
	sp->coefs->data[9] = 0.748151592823710;
	sp->coefs->data[10] = 0.450541598502498;
	sp->coefs->data[11] = 0.0838213779969326;
	sp->coefs->data[12] = 0.228976968716819;
	sp->coefs->data[13] = 0.913337361501670;
	sp->coefs->data[14] = 0.152378018969223;
	sp->coefs->data[15] = 0.825816977489547;
	sp->coefs->data[16] = 0.538342435260057;
	sp->coefs->data[17] = 0.996134716626886;
	sp->coefs->data[18] = 0.0781755287531837;
	sp->coefs->data[19] = 0.442678269775446;

	printf("Begin Call globaltrans  ***************************************************\n");
	globaltrans( nsp, sp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);
	printf("Call globaltrans  Over ***************************************************\n");
	//globaltrans(nsp, sp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);

	//globaltrans(nsp, sp, Vmax, Amax, Jmax, RVmax_A, RVmax_C, RAmax, RJmax);
	M_free(nsp->coefs);
	M_free(nsp->knots);
	free(nsp);
	M_free(sp->coefs);
	M_free(sp->knots);
	free(sp);



	// 结束测试 globaltrans


	/*
	int Count, * retWidth;
	long uu;
	char* rv;

	Count = 12;
	retWidth = 0;
	char* name;

	name = malloc(36);

	uu = getnew(Count, retWidth, name);
	
	printf("%ld\n", uu);
	printf("%s\n", name);
	system("pause");
	exit(0);
	

	// TODO: 在此添加控件通知处理程序代码
	MATRIX_TYPE _mat_1[3][5] = { 2, 2, 1, 1, 1, 4, 1, 1, 1, 1, 1, 5, 1, 1, 1 };
	int row = sizeof(_mat_1) / sizeof(_mat_1[0]);
	int column = sizeof(_mat_1[0]) / sizeof(_mat_1[0][0]);
	MATRIX_TYPE _mat_A10[10][10] = { 0.221746734017240, 0.928854139478045, 0.488897743920167, 0.0987122786555743,
									0.715037078400694, 0.904722238067363, 0.489901388512224, 0.521649842464284,
									0.800330575352402, 0.0604711791698936, 0.117417650855806, 0.730330862855453,
									0.624060088173690, 0.261871183870716, 0.903720560556316, 0.609866648422558,
									0.167927145682257, 0.0967300257808670, 0.453797708726920, 0.399257770613576,
									0.296675873218327, 0.488608973803579, 0.679135540865748, 0.335356839962797,
									0.890922504330789, 0.617666389588455, 0.978680649641159, 0.818148553859625,
									0.432391503783462, 0.526875830508296, 0.318778301925882, 0.578525061023439,
									0.395515215668593, 0.679727951377338, 0.334163052737496, 0.859442305646212,
									0.712694471678914, 0.817547092079286, 0.825313795402046, 0.416799467930787,
									0.424166759713807, 0.237283579771521, 0.367436648544477, 0.136553137355370,
									0.698745832334795, 0.805489424529686, 0.500471624154843, 0.722439592366842,
									0.0834698148589140, 0.656859890973707, 0.507858284661118, 0.458848828179931,
									0.987982003161633, 0.721227498581740, 0.197809826685929, 0.576721515614685,
									0.471088374541939, 0.149865442477967, 0.133171007607162, 0.627973359190104,
									0.0855157970900440, 0.963088539286913, 0.0377388662395521, 0.106761861607241,
									0.0305409463046367, 0.182922469414914, 0.0596188675796392, 0.659605252908307,
									0.173388613119006, 0.291984079961715, 0.262482234698333, 0.546805718738968,
									0.885168008202475, 0.653757348668560, 0.744074260367462, 0.239932010568717,
									0.681971904149063, 0.518594942510538, 0.390937802323736, 0.431651170248720,
									0.801014622769739, 0.521135830804002, 0.913286827639239, 0.494173936639270,
									0.500022435590201, 0.886511933076101, 0.0424311375007417, 0.972974554763863,
									0.831379742839070, 0.0154871256360190, 0.0292202775621463, 0.231594386708524,
									0.796183873585212, 0.779051723231275, 0.479922141146060, 0.0286741524641061,
									0.0714454646006424, 0.648991492712356, 0.803364391602440, 0.984063724379154 };
	row = sizeof(_mat_A10) / sizeof(_mat_A10[0]);
	column = sizeof(_mat_A10[0]) / sizeof(_mat_A10[0][0]);
	//Matrix* mat_A10 = Matrix_gen(row, column, _mat_A10);
	//M_print(mat_A10);



	Matrix *ncoef, *nt, * alpha, * beta;
	stusp1 sp1;

	//sp = (stusp1*)malloc(sizeof(stusp1*));

	sp1.number = 4;

	sp1.knots = Matrix_gen1(1, 8, 0.00);

	sp1.knots->row = 1;
	sp1.knots->column = 8;
	//sp1.knots->data = (MATRIX_TYPE*)malloc(sp1.knots->row * sp1.knots->column * sizeof(MATRIX_TYPE));
	sp1.knots->data[0] = 0;
	sp1.knots->data[1] = 0;
	sp1.knots->data[2] = 0;
	sp1.knots->data[3] = 0;
	sp1.knots->data[4] = 1;
	sp1.knots->data[5] = 1;
	sp1.knots->data[6] = 1;
	sp1.knots->data[7] = 1;


	sp1.coefs = Matrix_gen(sizeof(_mat_A10) / sizeof(_mat_A10[0]), sizeof(_mat_A10[0]) / sizeof(_mat_A10[0][0]), _mat_A10);
	//sp1.coefs->row = sizeof(_mat_A10) / sizeof(_mat_A10[0]);
	//sp1.coefs->column = sizeof(_mat_A10[0]) / sizeof(_mat_A10[0][0]);

	M_print( sp1.coefs);


	MATRIX_TYPE _nt1[1][12] = { 0.221746734017240, 0.928854139478045, 0.488897743920167, 0.0987122786555743,
								0.715037078400694, 0.904722238067363, 0.489901388512224, 0.521649842464284,
								0.800330575352402, 0.0604711791698936, 0.117417650855806, 0.730330862855453};


	nt = Matrix_gen(sizeof(_nt1) / sizeof(_nt1[0]), sizeof(_nt1[0]) / sizeof(_nt1[0][0]), _nt1);

	MATRIX_TYPE _alpha1[1][12] = { 0.33, 0.92, 0.488897743920167, 0.0987122786555743,
								0.715037078400694, 0.904722238067363, 0.489901388512224, 0.521649842464284,
								0.800330575352402, 0.0604711791698936, 0.117417650855806, 0.730330862855453 };

	//测试 spmak2
	alpha = Matrix_gen(sizeof(_alpha1) / sizeof(_alpha1[0]), sizeof(_alpha1[0]) / sizeof(_alpha1[0][0]), _alpha1);
	MATRIX_TYPE _beta1[1][12] = { 0.11, 0.928854139478045, 0.488897743920167, 0.0987122786555743,
								0.715037078400694, 0.904722238067363, 0.489901388512224, 0.521649842464284,
								0.800330575352402, 0.0604711791698936, 0.117417650855806, 0.730330862855453 };

	beta = Matrix_gen(sizeof(_beta1) / sizeof(_beta1[0]), sizeof(_beta1[0]) / sizeof(_beta1[0][0]), _beta1);


	//ncoef->row = 1;
	//ncoef->column = 1;

	MATRIX_TYPE* retdata;
	ncoef = Matrix_gen1(5, sp1.number, 0.00);

	stusp1* sp2;

	sp2 = (stusp1*)malloc(sizeof(stusp1*));

	sp2->number = 6;

	//M_print(sp2->knots);
	spmak2(alpha, beta, sp2);

	M_print(sp2->knots);

	//sp2->knots = beta;
	//sp2->coefs = alpha;

	//M_print(sp2->knots);

	//测试bangc
	bangc( ncoef, sp1, *nt, *alpha,  *beta);

	//测试  multal
	Matrix *tknot1;
	Matrix * nalpha1, * nbeta1;

	//nalpha = alpha * ones(1, n - 3);
	nalpha1 = Matrix_gen1(1, sp2->number - 3, 1 * 1.32);

	//nbeta = zeros(1, n - 3);
	nbeta1 = Matrix_gen1(1, sp2->number - 3, 0);

	//tknot = zeros(1, n + 4);
	tknot1 = Matrix_gen1(1, sp2->number - 4, 0);

	multal( tknot1, nalpha1, nbeta1, 1.32, sp2);

	//M_print(sp2->coefs);
	//M_print(sp2->knots);
	
	//M_print(tknot1);
	//M_print(nalpha1);
	//M_print(nbeta1);

	//ncoef->data = retdata;


	//printf("ncoef0= %lf\n", retdata[1]);

	//测试fnval2 

	stusp1 *sp3;

	sp3 = (stusp1*)malloc(sizeof(stusp1*));

	sp3->number = 4;

	sp3->knots = Matrix_gen1(1, 8, 0.00);

	sp3->knots->row = 1;
	sp3->knots->column = 3;
	//sp1.knots->data = (MATRIX_TYPE*)malloc(sp1.knots->row * sp1.knots->column * sizeof(MATRIX_TYPE));
	sp3->knots->data[0] = 0.11;
	sp3->knots->data[1] = 0.22;
	sp3->knots->data[2] = 0.33;
	*/
	/*
	sp3->knots->data[3] = 0.44;
	sp3->knots->data[4] = 1.11;
	sp3->knots->data[5] = 1.22;
	sp3->knots->data[6] = 1.33;
	sp3->knots->data[7] = 1.44;
	*/

  /*
	MATRIX_TYPE _mat_sp3c[3][4] = { 0.11, 0.22, 0.33, 0.44,
								0.55, 0.66, 0.77, 0.88,
								1.11, 1.22, 1.33, 1.44 };
	int row1, column1;

	row1 = sizeof(_mat_sp3c) / sizeof(_mat_sp3c[0]);
	column1 = sizeof(_mat_sp3c[0]) / sizeof(_mat_sp3c[0][0]);

	sp3->coefs = Matrix_gen(sizeof(_mat_sp3c) / sizeof(_mat_sp3c[0]), sizeof(_mat_sp3c[0]) / sizeof(_mat_sp3c[0][0]), _mat_sp3c);

	Matrix* v, *s;
	int m, r;
	m = 0;
	r = 0;
	m = sp3->coefs->row;


	MATRIX_TYPE _s[1][12] = { 0.11, 0.22, 0.33, 0.44,
								0.55, 0.66, 0.77, 0.88};

	s = Matrix_gen(sizeof(_beta1) / sizeof(_beta1[0]), sizeof(_beta1[0]) / sizeof(_beta1[0][0]), _s);

	

	r = s->column;

	v = Matrix_gen1(m, r, 0.00);

	fnval2a(v, sp3, s);



	//M_print(ncoef);

	M_free(tknot1);
	M_free(nalpha1);
	M_free(nbeta1);

	if (ncoef) 	
		M_free(ncoef);
	if (beta) 	
		M_free(beta);
	if (alpha) 	
		M_free(alpha);
	if (nt) 	
		M_free(nt);
	if (sp1.coefs) 	 
		M_free(sp1.coefs);
	if (sp1.knots->data) 	
		free(sp1.knots->data);
	*/

	//system("pause");
	
	//exit(0);






}
//
// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
