#include "datagroup.h"


/*
% ---- - 数据点分组
% ---- - 输出分组断点group，离散曲率k，数据点间距d
% ---- - 输入数据点data double型3 * n，曲率阈值kmax double型1 * 1，长直线长度dmax double型1 * 1
function[group, k, d] = datagroup(data, kmax, dmax)
*/

void datagroup(CellMat group, CellMat k, CellMat d, const CellMat data, const double kmax, const double dmax)
{

	//r = length(data(1, :));
	int r;

	r = 0;
	r = data.iCol;

	//% ---- - 计算离散曲率值---- -
	//d = ((data(1, 2:end) - data(1, 1:end - 1)). ^ 2 + (data(2, 2:end) - data(2, 1:end - 1)). ^ 2 + (data(3, 2:end) - data(3, 1:end - 1)). ^ 2). ^ (1 / 2); % || Pk + 1 - Pk ||
	//CellMat d;
	int i;
	i = 0;
	for (i=1; i<data.iCol-1; i++)
	{
		d.CellData[0][i - 1] = data.CellData[0][i] - data.CellData[0][i - 1];
		//data(1, 2:end) - data(1, 1:end - 1)
	}

	for (i = 0; i < data.iCol - 1; i++)
	{
		d.CellData[0][i] = d.CellData[0][i] * d.CellData[0][i];
	}

	for (i = 1; i < data.iCol - 1; i++)
	{
		d.CellData[1][i - 1] = data.CellData[1][i] - data.CellData[1][i - 1];
		//data(1, 2:end) - data(1, 1:end - 1)
	}

	for (i = 0; i < data.iCol - 1; i++)
	{
		d.CellData[1][i] = d.CellData[1][i] * d.CellData[1][i];
	}

	for (i = 1; i < data.iCol - 1; i++)
	{
		d.CellData[2][i - 1] = data.CellData[2][i] - data.CellData[2][i - 1];
		//data(1, 2:end) - data(1, 1:end - 1)
	}

	for (i = 0; i < data.iCol - 1; i++)
	{
		d.CellData[2][i] = d.CellData[2][i] * d.CellData[2][i];
	}

	for (i = 0; i < data.iCol - 1; i++)
	{
		d.CellData[0][i] = d.CellData[0][i] + d.CellData[1][i] + d.CellData[2][i];
		d.CellData[0][i] = sqrt( d.CellData[0][i] );
	}

	//d2 = ((data(1, 3:end) - data(1, 1:end - 2)). ^ 2 + (data(2, 3:end) - data(2, 1:end - 2)). ^ 2 + (data(3, 3:end) - data(3, 1:end - 2)). ^ 2). ^ (1 / 2);
	CellMat d2;

	for (i = 1; i < data.iCol - 2; i++)
	{
		d2.CellData[0][i - 1] = data.CellData[0][i+1] - data.CellData[0][i - 1];
	}

	for (i = 0; i < data.iCol - 2; i++)
	{
		d2.CellData[0][i] = d2.CellData[0][i] * d2.CellData[0][i];
	}

	for (i = 1; i < data.iCol - 2; i++)
	{
		d2.CellData[1][i - 1] = data.CellData[1][i] - data.CellData[1][i - 1];
	}

	for (i = 0; i < data.iCol - 2; i++)
	{
		d2.CellData[1][i] = d2.CellData[1][i] * d2.CellData[1][i];
	}

	for (i = 1; i < data.iCol - 2; i++)
	{
		d2.CellData[2][i - 1] = data.CellData[2][i] - data.CellData[2][i - 1];
	}

	for (i = 0; i < data.iCol - 2; i++)
	{
		d2.CellData[2][i] = d2.CellData[2][i] * d2.CellData[2][i];
	}

	for (i = 0; i < data.iCol - 2; i++)
	{
		d2.CellData[0][i] = d2.CellData[0][i] + d2.CellData[1][i] + d2.CellData[2][i];
		d2.CellData[0][i] = sqrt(d2.CellData[0][i]);
	}

	//p = (d(1:end - 1) + d(2:end) + d2) / 2;

	CellMat p;

	for (i = 0; i < data.iCol - 2; i++)
	{
		p.CellData[0][i] = (d.CellData[0][i] = d.CellData[1][i + 1] + d2.CellData[0][i]) / 2;
	}

	//s = (p.*(p - d(1:end - 1)).*(p - d(2:end)).*(p - d2)). ^ (1 / 2); % 面积
	//s = (p.*(p - d(1:end - 1)).*(p - d(2:end)).*(p - d2)). ^ (1 / 2);
	CellMat s, s1, s2, s3, s4, s5;

	for ( i=0; i < d.iCol-1; i++)
	{
		s1.CellData[0][i] = p.CellData[0][i] - d.CellData[0][i];
	}

	for (i = 1; i < d.iCol - 1; i++)
	{
		s2.CellData[0][i] = p.CellData[0][i] - d.CellData[0][i];
	}


	for (i = 0; i < p.iCol - 1; i++)
	{
		s3.CellData[0][i] = p.CellData[0][i] - d2.CellData[0][i];
	}

	for (i = 0; i < p.iCol - 1; i++)
	{
		//s3.CellData[0][i] = p.CellData[0][i] - d2.CellData[0][i];

		s.CellData[0][i] = sqrt(p.CellData[0][i] * s1.CellData[0][i] * s2.CellData[0][i] * s3.CellData[0][i]);

	}


	//k = [0, 2 * s.*(d(1:end - 1) + d(2:end)). ^ 2. / (d(1:end - 1).*d(2:end).*d2. ^ 3), 0];

	CellMat k0, kd1, kd2, kd3;

	for ( i=0; i < d.iCol-1; i++)
	{
		kd1.CellData[0][i] = d.CellData[0][i] + d.CellData[0][i+1];
	}

	for (i = 0; i < s.iCol - 1; i++)
	{
		kd2.CellData[0][i] = pow(2 * s.CellData[0][i] * kd1.CellData[0][i], 2); 
	}

	for (i = 0; i < d.iCol - 1; i++)
	{
		kd3.CellData[0][i] = pow(d.CellData[0][i] * d.CellData[0][i + 1] * d2.CellData[0][i], 3);
	}

	for (i = 0; i < kd3.iCol - 1; i++)
	{
		k0.CellData[0][i] = kd2.CellData[0][i] / kd3.CellData[0][i] ;
	}

	k.CellData[0][0] = 0;
	//k[1] = k0;
	for (i=0; i < k0.iCol; i++)
	{
		k.CellData[0][i + 1] = k0.CellData[0][i];

	}

	k.CellData[0][k0.iCol+1] = 0;


	//% ---- - 根据离散曲率分组---- -
	//group = zeros(r, 1);
	ZXMatInit(group, r, r, 1);

	//group(1) = 1;
	group.CellData[0][0] = 1;

	//group(r) = r;
	group.CellData[r%group.iRow][r/group.iRow] = r;

	/*
	for i = 1:r - 2 % 折线长度大于dmax的保留折线段
		if d(i) >= dmax % || d(i + 1) / d(i) >= 2 || d(i + 1) / d(i) <= 0.5
			group(i) = i;
			group(i + 1) = i + 1;
		end
	end
	*/

	for (i=1; i <= r-2; i++)
	{
		if (d.CellData[0][0] >= dmax)
		{
			group.CellData[i % group.iRow][i / group.iRow] = i;
			group.CellData[(i +1) % group.iRow][(i+1) / group.iRow] = i + 1;
		}
	}
	/*
	if d(r - 1) >= dmax
		group(r - 1) = r - 1;
		group(r) = r;
	end
	*/

	if (d.CellData[(r-1)%d.iRow][(r-1)/d.iRow] >= dmax)
	{
		group.CellData[(r - 1) % d.iRow][(r - 1) / d.iRow] = r - 1;
		group.CellData[r % d.iRow][r / d.iRow] = r ;
	}

	/*
	for i = 2:r - 1 % 根据曲率增加分段点
		k_ = abs(k);
		if k_(i) >= kmax % || k_(i + 1) / k_(i) >= 4 || k_(i + 1) / k_(i) <= 0.25
			group(i) = i;
		end
	end
	*/

	for (i=2; i <= r-1; i++)
	{

		if (k.CellData[i%k.iRow][i/k.iRow] >= kmax )
		{
			group.CellData[i % k.iRow][i / k.iRow] = i;
		}


	}

	//	% group(group == 0) = [];% 分段点
	//group = [group(1:end - 1)';group(2:end)'];

	CellMat  groupTmp;

	groupTmp.iRow = 2;
	groupTmp.iCol = group.iCol - 1;

	for (i=0; i< groupTmp.iCol; i++)
	{
		groupTmp.CellData[0][i] = group.CellData[0][i];
		groupTmp.CellData[1][i] = group.CellData[0][i+1];
	}


	group.iRow = groupTmp.iRow;
	group.iCol = groupTmp.iCol;

	for (i = 0; i < groupTmp.iCol; i++)
	{
		group.CellData[0][i] = groupTmp.CellData[0][i];
		group.CellData[1][i] = groupTmp.CellData[0][i];
	}

}


