function [E,L] = eigwithknown0(M)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
MT = M';
v3 = crossProduct(MT(1,:),MT(2,:));
v3 = normalizeVector(v3);
% L(1,3) = 0;
% v3 = [x(2,1)*x(3,2)-x(3,1)*x(2,2),x(3,1)*x(1,2)-x(3,2)*x(1,1),x(2,2)*x(1,1)-x(2,1)*x(1,2)];
% v3 = v3 / norm(v3);

% x01_squared = SQR(M(1,2));
b = -M(1,1) - M(2,2) - M(3,3);
c = -SQR(M(1,2)) - SQR(M(1,3)) - SQR(M(2,3)) + M(1,1)*(M(2,2) + M(3,3)) + M(2,2)*M(3,3);
[e1,e2] = root2real(b,c);
if (abs(e1) < abs(e2))
    temp = e1;
    e1 = e2;
    e2 = temp;
end
L(1,1) = e1;
L(1,2) = e2;
L(1,3) = 0;

mx0011 = -M(1,1)*M(2,2);
prec_0 = M(1,2)*M(2,3) - M(1,3)*M(2,2);
prec_1 = M(1,2)*M(1,3) - M(1,1)*M(2,3);

tmp = 1 / (e1*(M(1,1) + M(2,2)) + mx0011 - SQR(e1) + SQR(M(1,2)));
a1 = -(e1*M(1,3) + prec_0)*tmp;
a2 = -(e1*M(2,3) + prec_1)*tmp;
rnorm = 1 / sqrt(SQR(a1) + SQR(a2) + 1);
a1 = a1*rnorm;
a2 = a2*rnorm;
v1 = [a1,a2,rnorm];

tmp2 = 1 / (e2*(M(1,1) + M(2,2)) + mx0011 - SQR(e2) + SQR(M(1,2)));
a21 = -(e2*M(1,3) + prec_0)*tmp2;
a22 = -(e2*M(2,3) + prec_1)*tmp2;
rnorm2 = 1 / sqrt(SQR(a21) + SQR(a22) + 1);
a21 = a21 * rnorm2;
a22 = a22 * rnorm2;
v2 = [a21,a22,rnorm2];

E = [v1(1,1),v2(1,1),v3(1,1);v1(1,2),v2(1,2),v3(1,2);v1(1,3),v2(1,3),v3(1,3)];


end

