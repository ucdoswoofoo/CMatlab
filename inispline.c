#include "inispline.h"
#include "featurepoint.h"
#include "spmak2.h"


/*
% �������ݵ㣬�����ʼ������ߺ����
function sp = inispline(segdata, segwrefq, eps)

m = length(segdata(1, :));

% ��������ϵ��������
u = zeros(1, m);
for i = 2:m
u(i) = u(i - 1) + norm(segwrefq(:, i) - segwrefq(:, i - 1));
end
u = u / u(m);
*/

void inispline(stusp1* sp, Matrix segdata, Matrix segwrefq, const double eps)
{
	//m = length(segdata(1, :));

	
}


