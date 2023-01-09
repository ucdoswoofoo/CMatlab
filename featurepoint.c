#include "featurepoint.h"
#include "spmak2.h"
#include "data3curve.h"

/*
% ---- - 选择特征点
% ---- - 输出特征点位置（索引）feature，离散曲率k，数据点间距d
% ---- - 输入数据点data，拟合精度dmax
function[feature, k, d] = featurepoint(segwrefq, dmax)
r = length(segwrefq(1, :));
% ---- - 计算离散曲率值---- -

*/
	// featurepoint(CellMat feature, double k, CellMat d, CellMat* segwrefq, const double  dmax);
void featurepoint(Matrix feature, Matrix k, Matrix d, Matrix* segwrefq, const double dmax)
{
	int r;
	r = 0;
	
}


