%Author: Yu, Qida
%This programe generates experiment data;

function [XXw,xn,R,t]=genP3Pdata(noise)

npt=3;
% camera's parameters
width= 640;
height= 480;
f=800;
nl=noise;
% generate 3d coordinates in camera space
Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
Xcc= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)];
R = getRotationMatrix(getRandomUnitVector(4));
t= mean(Xc,2);
XXw= (R)\(Xc-repmat(t,1,npt));
xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
xxn= xx+randn(2,npt)*nl;
%normalized focal;
xn=xxn./f;
iscolinear = 0;
for i = 1:3
    for j = 1:3
        if(i ~= j)
            if (((XXw(1,i)==XXw(1,j))&&(XXw(2,i)==XXw(2,j))&&(XXw(3,i)==XXw(3,j)))|| ((Xcc(1,i)==Xcc(1,j))&&(Xcc(2,i)==Xcc(2,j))))
                iscolinear = 1;
            end
        end
    end
end
if (iscolinear)
     [XXw,xn,R,t]=genP3Pdata(noise);
end
end
