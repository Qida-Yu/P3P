%Author: Yu, Qida
%This is our codes to solve P3P problem.

function res = P3P_Yu(P,xn)

m1=xn(:,1); 
m2=xn(:,2); 
m3=xn(:,3); 
P1=P(:,1);
P2=P(:,2);
P3=P(:,3);
ux=(P2-P1)/norm(P2-P1);
uz=cross(ux,P3-P1)/norm(cross(ux,P3-P1));
uy=cross(uz,ux);
N=[ux,uy,uz]';
tM = [zeros(3,1),P2-P1,P3-P1];
tM = N * tM;
a = tM(1,2);
b = tM(1,3);
c = tM(2,3);

b10 = (m2(1) - m1(1)) / a;
b30 = (m2(2) - m1(2)) / a;

b21 = (m3(1) - m2(1)) * b / c;
b20 = (m3(1) - m1(1) - b * b10) / c;
b41 = (m3(2) - m2(2)) * b / c;
b40 = (m3(2) - m1(2) - b * b30) / c;

g0=b21*m2(1)+b41*m2(2);
g1=m2'*m3+1;
g2=m2(1)*b20+b21*b10+b41*b30+m2(2)*b40;
g3=m3(1)*b10+m3(2)*b30;
g4=b20*b10+b40*b30;

h0 = SQR(b21) + SQR(b41) - m2'*m2 - 1;
h1 = 2*m3(1)*b21+2*m3(2)*b41;
h2 = SQR(m3(1)) + SQR(m3(2)) + 1;
h3 = 2*b21*b20+2*b41*b40-2*m2(1)*b10-2*m2(2)*b30;
h4 = 2*m3(1)*b20+2*m3(2)*b40;
h5 = SQR(b20) + SQR(b40) - SQR(b10) - SQR(b30);

k0=h2*SQR(g0)+h0*SQR(g1)-h1*g0*g1;
k1=2*h2*g0*g2+2*h0*g1*g3-h1*(g0*g3+g1*g2)-h4*g0*g1+h3*SQR(g1);
k2=h2*SQR(g2)+2*h2*g0*g4+h0*SQR(g3)-h1*(g2*g3+g1*g4)-h4*(g0*g3+g1*g2)+2*h3*g1*g3+h5*SQR(g1);
k3=2*h2*g2*g4-h1*g3*g4-h4*(g2*g3+g1*g4)+h3*SQR(g3)+2*h5*g1*g3;
k4=h2*SQR(g4)-h4*g3*g4+h5*SQR(g3);
U=[k0,k1,k2,k3,k4];

a_roots  = roots(U);
a_roots1 = polishRoots(U,a_roots);
a_roots1 = real(a_roots);
rindex = abs(a_roots1) <= 1/4;
a_roots1 = a_roots1(rindex);
n1=length(a_roots1);

for i = 1 : n1
    alpha1 = a_roots1(i);
    alpha2 = -(g0*SQR(alpha1)+g2*alpha1+g4)/(g1*alpha1+g3);

    nrm = SQR(alpha1*m2(1)+b10)+SQR(alpha1*m2(2)+b30)+SQR(alpha1)+SQR(alpha1*b21+alpha2*m3(1)+b20)+SQR(alpha1*b41+alpha2*m3(2)+b40)+SQR(alpha2);
    s = 1 / sqrt(nrm/2);
    r1 = (alpha1*m2(1)+b10)*s;
    r2 = (alpha1*b21+alpha2*m3(1)+b20)*s;
    r4 = (alpha1*m2(2)+b30)*s;
    r5 = (alpha1*b41+alpha2*m3(2)+b40)*s;
    r7 = alpha1*s;
    r8 = alpha2*s;
    r3 = r4 * r8 - r5 * r7;
    r6 = r2 * r7 - r1 * r8;
    r9 = r1 * r5 - r2 * r4;
    tx = m1(1)*s;
    ty = m1(2)*s;
    T_F=[tx;ty;s];
    R_F=[r1,r2,r3;r4,r5,r6;r7,r8,r9];
    res{i}.R=R_F*N;
    res{i}.T=T_F-R_F*N*P1;
end

end

