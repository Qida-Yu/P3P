%Author: Yu, Qida
%This is the codes to test the computational time.

close all;
clc;
clear;

addpath solutions;
addpath functions;

N=200;
noise=0;

%generate data;
for i=1:N
    [Xw,xn,R,t]=genP3Pdata(noise);
    t1=clock;
    for j=1:500 
         res1=P3P_Gao(Xw,xn);
     end
     t2=clock;
     time1(i)=etime(t2,t1);
     
     t3=clock;
     for j=1:500  
         res2=P3P_K(Xw,xn);
    end
     t4=clock;
     time2(i)=etime(t4,t3);
     
     t5=clock;
     for j=1:500 
         res3=P3P_B(Xw,xn);
     end
     t6=clock;
     time3(i)=etime(t6,t5);
     
     t7=clock;
     for j=1:500 
         res4=P3P_Ke(Xw,xn);
     end
     t8=clock;
     time4(i)=etime(t8,t7);
     
     t9=clock;
     for j=1:500 
         res5=P3P_Persson(Xw,xn);
     end
     t10=clock;
     time5(i)=etime(t10,t9);
    
    t11=clock;
    for j=1:500 
        res6=P3P_Yu(Xw,xn);
    end
    t12=clock;
    time6(i)=etime(t12,t11);

    t13=clock;
    for j=1:500 
        res7=P3P_Jrrrr(Xw,xn);
    end
    t14=clock;
    time7(i)=etime(t14,t13);

 sum(time1)
 sum(time2)
 sum(time3)
 sum(time4)
 sum(time5)
 sum(time6)
 sum(time7)

 Gao = sum(time1)/(N*500);
 Kneip = sum(time2)/(N*500);
 Banno = sum(time3)/(N*500);
 Ke = sum(time4)/(N*500);
 Persson = sum(time5)/(N*500);
 Yu= sum(time6)/(N*500);
 Jrrrr= sum(time7)/(N*500);
end