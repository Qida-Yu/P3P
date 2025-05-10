% codes for original paper
function res = P3P_PerssonNP(P,xn)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
x=[xn;ones(1,3)];
x = xnorm(x);
y1 = x(:,1)';
y2 = x(:,2)';
y3 = x(:,3)'; 
P = P';
x1 = P(1,:);
x2 = P(2,:);
x3 = P(3,:);

b12 = -2 * (y1*y2');
b13 = -2 * (y1*y3');
b23 = -2 * (y2*y3');

d12 = x1 - x2;
d13 = x1 - x3;
d23 = x2 - x3;
d12xd13 = cross(d12,d13);

a12 = d12*d12';
a13 = d13*d13';
a23 = d23*d23';

c31 = -0.5 * b13;
c23 = -0.5 * b23;
c12 = -0.5 * b12;
blob = (c12*c23*c31 - 1);

s31_squared = 1 - SQR(c31);
s23_squared = 1 - SQR(c23);
s12_squared = 1 - SQR(c12);

p3 = (a13*(a23*s31_squared - a13*s23_squared));
p2 = 2*blob*a23*a13 + a13*(2*a12 + a13)*s23_squared + a23*(a23 - a12)*s31_squared;
p1 = a23*(a13 - a23)*s12_squared - a12*a12*s23_squared - 2*a12*(blob*a23 + a13*s23_squared);
p0 = a12*(a12*s23_squared - a23*s12_squared);

g = cubick(p2/p3,p1/p3,p0/p3);

A00 = a23*(1 - g);
A01 = (a23*b12)*0.5;
A02 = (a23*b13*g)*(-0.5);
A11 = a23 - a12 + a13*g;
A12 = b23*(a13*g - a12)*(0.5);
A22 = g*(a13 - a23) - a12;

A = [A00,A01,A02;A01,A11,A12;A02,A12,A22];
[V,L] = eigwithknown0(A);
v = sqrt(max(0,-L(1,2)/L(1,1)));
valid1 = 0;


s = v;
w2 = 1/(s*V(1,2) - V(1,1));
w0 = (V(2,1) - s*V(2,2)) * w2;
w1 = (V(3,1) - s*V(3,2)) * w2;

a = 1/((a13 - a12)*SQR(w1) - a12*b13*w1 - a12);
b = (a13*b12*w1 - a12*b13*w0 - 2*w0*w1*(a12 - a13)) * a;
c = ((a13 - a12)*SQR(w0) + a13*b12*w0 + a13) * a;
tmp = SQR(b) - 4*c;
if (tmp >= 0)
    [tau1,tau2] = root2real(b,c);
    if (tau1 > 0)
        tau = tau1;
        d = a23 / (tau*(b23 + tau) + 1);
        l2 = sqrt(d);
        l3 = tau*l2;
        l1 = w0*l2 + w1*l3;
        if (l1 >= 0)
            valid1 = valid1 + 1;
            Ls(valid1,:) = [l1,l2,l3];
        end
    end
    if (tau2 > 0)
        tau = tau2;
        d = a23 / (tau*(b23 + tau) + 1);
        l2 = sqrt(d);
        l3 = tau*l2;
        l1 = w0*l2 + w1*l3;
        if (l1 >= 0)
            valid1 = valid1 + 1;
            Ls(valid1,:) = [l1,l2,l3];
        end
    end
end

s = -v;
w2 = 1 / (s*V(1,2) - V(1,1));
w0 = (V(2,1) - s*V(2,2)) * w2;
w1 = (V(3,1) - s*V(3,2)) * w2;

a = 1/((a13 - a12)*SQR(w1) - a12*b13*w1 - a12);
b = (a13*b12*w1 - a12*b13*w0 - 2*w0*w1*(a12 - a13)) * a;
c = ((a13 - a12)*SQR(w0) + a13*b12*w0 + a13) * a;
tmp = (SQR(b) - 4*c);
if (tmp >= 0)
    [tau1,tau2] = root2real(b,c);
    if (tau1 > 0)
        tau = tau1;
        d = a23 / (tau*(b23 + tau) + 1);
          l2 = sqrt(d);
          l3 = tau*l2;
          l1 = w0*l2 + w1*l3;
          if (l1 >= 0)
             valid1 = valid1 + 1;
             Ls(valid1,:) = [l1,l2,l3];
          end
    end
    if (tau2 > 0)
        tau = tau2;
        d = a23 / (tau*(b23 + tau) + 1);
          l2 = sqrt(d);
          l3 = tau*l2;
          l1 = w0*l2 + w1*l3;
          if (l1 >= 0)
             valid1 = valid1 + 1;
             Ls(valid1,:) = [l1,l2,l3];
          end
    end
end

% for i = 1 : valid1
%     Ls(i,:) = gauss_newton_refineL(Ls(i,:),a12,a13,a23,b12,b13,b23);
% end

X = [d12',d13',d12xd13'];

for j = 1 : valid1
    ry1 = y1*Ls(j,1);
    ry2 = y2*Ls(j,2);
    ry3 = y3*Ls(j,3);
    
    yd1 = ry1 - ry2;
    yd2 = ry1 - ry3;
    yd1xd2 = cross(yd1,yd2);
    Y = [yd1',yd2',yd1xd2'];
    R_tmp =  Y / X;
    res{j}.Xc(:,1)=ry1';
    res{j}.Xc(:,2)=ry2';
    res{j}.Xc(:,3)=ry3';
    res{j}.R= R_tmp;
    res{j}.T= ry1' - R_tmp*x1';
end
    
end

