function root = cubick_new(c2,c1,c0)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
a = c1 - SQR(c2) / 3;
b = (2 * TRI(c2) - 9 * c1 * c2) / 27 + c0;
c = SQR(b) / 4 + TRI(a) / 27;
if (c > 1.0e-26) 
    c = sqrt(c);
    b = -0.5*b;
    root = nthroot((b + c),3) + nthroot((b - c),3) - c2 / 3;
elseif (c < -1.0e-26) 
    c = 3 * b / (2 * a) * sqrt(-3 / a);
    root = 2 * sqrt(-a / 3.0) * cos(acos(c) / 3) - c2 / 3;
elseif (c <= 1.0e-26 && c >= -1.0e-26 && a~=0)
    root = 3*b/a - c2 / 3;
else 
    root = -c2 / 3;
end



