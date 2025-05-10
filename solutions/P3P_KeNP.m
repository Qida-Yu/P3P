%Author: Yu, Qida
%This programe to solve the P3P problem as described in: 
%T. Ke and S. I. Roumeliotis, ¡°An efficient algebraic
%solution to the perspective-three-point problem,¡± in
%2017 IEEE Conference on Computer Vision and Pattern
%Recognition (CVPR), 2017, pp. 4618¨C4626.
%Ke's method without root polishing

function res = P3P_KeNP(P,xn)

x=[xn;ones(1,3)];
b1=x(:,1)'./norm(x(:,1)); 
b2=x(:,2)'./norm(x(:,2)); 
b3=x(:,3)'./norm(x(:,3)); 
P1=P(:,1)';
P2=P(:,2)';
P3=P(:,3)';
u0 = P1 - P2;
nu0 = norm(u0);
k1 = u0 / nu0;
k3 = crossProduct(b1,b2);
nk3 = norm(k3);
k3 = k3 / nk3;
tz = crossProduct(b1,k3);
v1 = crossProduct(b1,b3);
v2 = crossProduct(b2,b3);

u1 = P1 - P3;
u1k1 = DOT(u1,k1);
k3b3 = DOT(k3,b3);

f11 = k3b3;
f13 = DOT(k3,v1);
f15 = -u1k1 * f11;

nl1 = crossProduct(u1,k1);
delta = norm(nl1);
nl1 = nl1 / delta;
f11 = f11 * delta;
f13 = f13 * delta;

u2k1 = u1k1 - nu0;
f21 = DOT(tz,v2);
f22 = nk3 * k3b3;
f23 = DOT(k3,v2);
f24 = u2k1 * f22;
f25 = -u2k1 * f21;
f21 = f21 * delta;
f22 = f22 * delta;
f23 = f23 * delta;

g1 = f13 * f22;
g2 = f13 * f25 - f15 * f23;
g3 = f11 * f23 - f13 * f21;
g4 = -f13 * f24;
g5 = f11 * f22;
g6 = f11 * f25 - f15 * f21;
g7 = -f15 * f24;

factor_4 = SQR(g5) + SQR(g1) + SQR(g3);
factor_3 = 2 * (g5 * g6 + g1 * g2 + g3 * g4);
factor_2 = SQR(g6) + 2 * g5 * g7 + SQR(g2) + SQR(g4) - SQR(g1) - SQR(g3);
factor_1 = 2 * (g6 * g7 - g1 * g2 - g3 * g4);
factor_0 = SQR(g7) - SQR(g2) - SQR(g4);

U = [factor_4, factor_3, factor_2, factor_1, factor_0];

% s_root = SolveQuartic(U);
s_root = roots(U);
s_roots = real(s_root);
rindex = abs(s_roots) <= 1;
s_roots1 = s_roots(rindex);
temp = crossProduct(k1,nl1);
Ck1nl = [k1;nl1;temp]';
Cb1k3tzT = [b1;k3;tz];
b3p = (delta / k3b3) * b3;
len1 = length(s_roots1);

for i = 1 : len1
    ctheta1p = s_roots1(i);
    stheta1p = sqrt(1 - SQR(ctheta1p));
    if (k3b3 <= 0)
        stheta1p = -stheta1p;
    end
    ctheta3 = g1 * ctheta1p + g2;
    stheta3 = g3 * ctheta1p + g4;
    ntheta3 = stheta1p / ((g5 * ctheta1p + g6) * ctheta1p + g7);
    ctheta3 = ctheta3 * ntheta3;
    stheta3 = stheta3 * ntheta3;
    
    C13 = [ctheta3, 0, -stheta3; stheta1p * stheta3, ctheta1p, stheta1p * ctheta3; ctheta1p * stheta3, -stheta1p, ctheta1p * ctheta3];
    
    temp_matrix = Ck1nl * C13;
    R1 = temp_matrix * Cb1k3tzT; 
    rp3 = [P3(1,1)*R1(1,1)+P3(1,2)*R1(2,1)+P3(1,3)*R1(3,1),P3(1,1)*R1(1,2)+P3(1,2)*R1(2,2)+P3(1,3)*R1(3,2),P3(1,1)*R1(1,3)+P3(1,2)*R1(2,3)+P3(1,3)*R1(3,3)];
    pxstheta1p = stheta1p * b3p;
%     the following is not stable for computing t
%     p = R1 * b3p'; 
%     pxstheta1p = stheta1p * p;
%     temp = P3' - pxstheta1p;
    res{i}.R = R1';
    res{i}.T = (pxstheta1p - rp3)';    
end

end

