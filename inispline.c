#include "inispline.h"
#include "featurepoint.h"
#include "spmak2.h"


/*
% 输入数据点，输出初始拟合曲线和误差
function sp = inispline(segdata, segwrefq, eps)

m = length(segdata(1, :));

% 工件坐标系弧长参数
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


