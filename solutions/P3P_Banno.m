%Author: Yu, Qida
%This programe to solve the P3P problem as described in: 
%A. Banno, ¡°A p3p problem solver representing all parameters
%as a linearcombination,¡± Image and Vision
%Computing, vol. 70, pp. 55 ¨C 62, 2018.


function res = P3P_Banno(P,xn)
x=[xn;ones(1,3)];
f1=x(:,1)/norm(x(:,1)); 
f2=x(:,2)/norm(x(:,2)); 
f3=x(:,3)/norm(x(:,3)); 
P1=P(:,1);
P2=P(:,2);
P3=P(:,3);
ux=(P2-P1)/norm(P2-P1);
uz=cross(ux,P3-P1)/norm(cross(ux,P3-P1));
uy=cross(uz,ux);
R_temp=[ux,uy,uz]';
P_u1=[0;0;0];
P_u2=R_temp*(P2-P1);
P_u3=R_temp*(P3-P1);
x_2=P_u2(1);
x_3=P_u3(1);
y_3=P_u3(2);
%----------------------
b_u1=zeros(12,1);
b_u2=zeros(12,1);
b_u3=zeros(12,1);
b_u1(10) = 1;
b_u2(11) = 1;
b_u3(12) = 1;

A1 = 1/x_2/y_3;
B1 = (x_2-x_3) * A1;

b_u1(1)= -f1(1) / x_2;
b_u1(2)= -f1(1) * B1;
b_u1(3)= -f1(2) / x_2;
b_u1(4)= -f1(2) * B1;
b_u1(5)= -f1(3) / x_2;
b_u1(6)= -f1(3) * B1;
b_u1(7)= f1(1);
b_u1(8)= f1(2);
b_u1(9)= f1(3);

b_u2(1)= f2(1) / x_2;
b_u2(2)= -f2(1) * x_3 * A1;
b_u2(3)= f2(2) / x_2;
b_u2(4)= -f2(2) * x_3 * A1;
b_u2(5)= f2(3) / x_2;
b_u2(6)= -f2(3) * x_3 * A1;

b_u3(2)= f3(1) / y_3;
b_u3(4)= f3(2) / y_3;
b_u3(6)= f3(3) / y_3;

c12 = f1'* f2;
c23 = f2'* f3;
c31 = f3'* f1;

y_32 = SQR(y_3);
r_12 = x_2 - x_3;

h1 = y_32 - r_12 * r_12;
h2 = -2.0 * (y_32 + r_12 * x_3) * c12;
h3 = 2.0 * x_2 * r_12 * c31;
h4 = y_32 - x_3 * x_3;
h5 = 2.0 * x_2 * x_3 * c23;
h6 = -x_2 * x_2;

g1 = r_12;
g2 = (x_3 - r_12) * c12;
g3 = -x_2 * c31;
g4 = -x_3;
g5 = x_2 * c23;

f11 = h1 * h1;
f12 = h1 * h2;
f13 = h1 * h3;
f14 = h1 * h4;
f15 = h1 * h5;
f16 = h1 * h6;

f22 = h2 * h2;
f23 = h2 * h3;
f24 = h2 * h4;
f25 = h2 * h5;
f26 = h2 * h6;

f33 = h3 * h3;
f34 = h3 * h4;
f35 = h3 * h5;
f36 = h3 * h6;

f44 = h4 * h4;
f45 = h4 * h5;
f46 = h4 * h6;

f55 = h5 * h5;
f56 = h5 * h6;

f66 = h6 * h6;

g11 = g1 * g1;
g12 = g1 * g2;
g13 = g1 * g3;
g14 = g1 * g4;
g15 = g1 * g5;

g22 = g2 * g2;
g23 = g2 * g3;
g24 = g2 * g4;
g25 = g2 * g5;

g33 = g3 * g3;
g34 = g3 * g4;
g35 = g3 * g5;

g44 = g4 * g4;
g45 = g4 * g5;

g55 = g5 * g5;

factor_4 = f44*g11 + f11*g44 - f24*g12 + f22*g14 + f14*g22 - f12*g24 - 2*f14*g14;
factor_3 = 2*(f45*g11 + f11*g45 + f23*g14 + f14*g23 - f15*g14 - f14*g15) - f34*g12 - f25*g12 - f24*g13 + f22*g15 + f15*g22 - f13*g24 - f12*g34 - f12*g25;
factor_2 = 2*(f46*g11 + f23*g15 + f15*g23 - f16*g14 - f15*g15) + f55*g11 + f11*g55 - f35*g12 - f26*g12 - f34*g13 - f25*g13 + f33*g14 + f16*g22 + f14*g33 - f13*g34 - f13*g25 - f12*g35;
factor_1 = 2*(f56*g11 + f16*g23 - f16*g15) - f36*g12 - f35*g13 - f26*g13 + f33*g15 + f15*g33 - f13*g35;
factor_0 = f66*g11 - f36*g13 + f16*g33;

sol_roots = roots([factor_4, factor_3, factor_2, factor_1, factor_0]);
sol_roots = real(sol_roots);
len1 = length(sol_roots);
for i = 1:len1
    solb = sol_roots(i);
    sola = h6 * g1 + solb * (h5 * g1 - h1 * g5 + solb * (h4 * g1 - h1 * g4));
    sola = sola / ((h1 * g2 - h2 * g1) * solb + h1 * g3 - h3 * g1);
    nrm = SQR(sola * b_u1(1) + solb * b_u2(1)) + SQR(sola * b_u1(3) + solb * b_u2(3)) + SQR(sola * b_u1(5) + solb * b_u2(5))...
        + SQR(sola * b_u1(2) + solb * b_u2(2) + b_u3(2)) + SQR(sola * b_u1(4) + solb * b_u2(4) + b_u3(4)) + SQR(sola * b_u1(6) + solb * b_u2(6) + b_u3(6));
    solc = 1/ sqrt(nrm/2);
    sola = sola * solc;
    solb = solb * solc;
    r1 = sola * b_u1(1) + solb * b_u2(1);
    r2 = sola * b_u1(2) + solb * b_u2(2) + solc * b_u3(2);
    r4 = sola * b_u1(3) + solb * b_u2(3);
    r5 = sola * b_u1(4) + solb * b_u2(4) + solc * b_u3(4);
    r7 = sola * b_u1(5) + solb * b_u2(5);
    r8 = sola * b_u1(6) + solb * b_u2(6) + solc * b_u3(6);
    r3 = r4 * r8 - r5 * r7;
    r6 = r2 * r7 - r1 * r8;
    r9 = r1 * r5 - r2 * r4;
    R_F = [r1,r2,r3;r4,r5,r6;r7,r8,r9];
    T_F=[sola * b_u1(7);sola * b_u1(8);sola * b_u1(9)];
    res{i}.R=R_F*R_temp;
    res{i}.T=T_F-R_F*R_temp*P1;
end
   
end

