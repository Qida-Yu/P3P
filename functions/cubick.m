function g = cubick(b,c,d)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
if(SQR(b) >= 3*c)
    v = sqrt(SQR(b) - 3*c);
    t1 = (-b - v)/3;
    
    k = ((t1 + b)*t1 + c)*t1 + d;
    if (k > 0)
        r0 = t1 - sqrt(-k / (3*t1 + b));
    else
        t2 = (-b + v) / 3;
        k = ((t2 + b)*t2 + c)*t2 + d;
        r0 = t2 + sqrt(-k / (3*t2 + b));
    end
else
    r0= -b / 3;
    if (abs(((3*r0+2*b)*r0+c)) < 1e-4)
        r0 = r0 + 1;
    end
end

for cnt = 1 : 50
    fx=(((r0 + b)*r0 + c)*r0+d);
    if (((cnt < 7) || (abs(fx) > 1e-15)))
        fpx = (3*r0 + 2*b)*r0 + c;
        r0 = r0 - fx / fpx;
    else
        break
    end
end
g = r0;
end

