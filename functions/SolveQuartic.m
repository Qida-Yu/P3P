%Author: Yu, Qida

function r = SolveQuartic(c)
% Ferrari's method to solve a quartic equation
a4 = c(1,1);
a3 = c(1,2);
a2 = c(1,3);
a1 = c(1,4);
a0 = c(1,5);

a4_2 = SQR(a4);
a3_2 = SQR(a3);
a4_3 = a4_2 * a4;
a2a4 = a2 * a4;

p4 = (8 * a2a4 - 3 * a3_2) / (8 * a4_2);
q4 = (a3_2 * a3 - 4 * a2a4 * a3 + 8 * a1 * a4_2) / (8 * a4_3);
r4 = (256 * a0 * a4_3 - 3 * (a3_2 * a3_2) - 64 * a1 * a3 * a4_2 + 16 * a2a4 * a3_2) / (256 * (a4_3 * a4));

p3 = (SQR(p4) / 12 + r4) / 3;
q3 = (72 * r4 * p4 - 2 * p4 * p4 * p4 - 27 * SQR(q4)) / 432;

if (q3 >= 0)
    w = -sqrt(q3 * q3 - p3 * p3 * p3) - q3;
else
    w = sqrt(q3 * q3 - p3 * p3 * p3) - q3;
end

if (imag(w) == 0.0)
    w_crbt = real(w);
    w_li = w_crbt^(1/3);
    t = 2.0 * (w_li + p3 / w_li);
else
    w = w^(1.0/3);
    t = 4.0 * real(w);
end

sqrt_2m = sqrt(-2 * p4 / 3 + t);
B_4A = -a3 / (4 * a4);
complex1 = 4 * p4 / 3 + t;
complex2 = 2 * q4 / sqrt_2m;
sqrt_2m_rh = real(sqrt_2m) / 2;
sqrt1 = real(sqrt(-(complex1 + complex2))) / 2;
r(1) = B_4A + sqrt_2m_rh + sqrt1;
r(2) = B_4A + sqrt_2m_rh - sqrt1;
sqrt2 = real(sqrt(-(complex1 - complex2))) / 2;
r(3) = B_4A - sqrt_2m_rh + sqrt2;
r(4) = B_4A - sqrt_2m_rh - sqrt2;
r = r';

end

