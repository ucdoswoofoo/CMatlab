#include "featurepoint.h"
#include "spmak2.h"
#include "data3curve.h"

/*
% ---- - ѡ��������
% ---- - ���������λ�ã�������feature����ɢ����k�����ݵ���d
% ---- - �������ݵ�data����Ͼ���dmax
function[feature, k, d] = featurepoint(segwrefq, dmax)
r = length(segwrefq(1, :));
% ---- - ������ɢ����ֵ---- -

*/
	// featurepoint(CellMat feature, double k, CellMat d, CellMat* segwrefq, const double  dmax);
void featurepoint(Matrix feature, Matrix k, Matrix d, Matrix* segwrefq, const double dmax)
{
	int r;
	r = 0;
	
}


