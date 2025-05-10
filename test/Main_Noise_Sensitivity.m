%This is the codes to test noise sensitivity.
close all;
clc;
clear;
addpath solutions;
addpath functions;


%The repetitive numbers to per noise level;
iterations=10000;
%The maximum noise to analyze
max_noise=5;
%The step in between dirrerent noise levels;
noise_step=0.1;
number_noise_levels=max_noise/noise_step + 1;
noise_levels=zeros(1,number_noise_levels);
%------------------------Algorithms comparison--------------------
%prepare the overall result arrays
C= zeros(1,iterations);
D= zeros(1,number_noise_levels);
name= {'Gao','Kneip','Banno','Ke','Persson','Yu', 'Ours'};
f= {@P3P_Gao,@P3P_Kneip,@P3P_Banno,@P3P_Ke,@P3P_Persson,@P3P_Yu,@P3P_Jrrr};
color = {'r',[0,0.7,0],'b',[1,0.5,0.1],'c','m',[0.5,0.2,0.8]}; 

method_list = struct('name', name, 'f', f, 'color', color, 'r1', C, 't1', C, 'mean_r1', D, 'mean_t1', D);



%beging;
for n=1:number_noise_levels
    noise=(n-1)*noise_step;
    noise_levels(n)=noise;
    display(['Analyzing noise level: ',num2str(noise)]);   
    
    for i=1:iterations
        %generating experiment data;
        [Xw,xn,R,t]=genP3Pdata(noise);
        for k=1:length(method_list)
            try
                res= method_list(k).f(Xw,xn);
            catch
                continue;
            end
            %Calculation errors;
            [res_out1,position_error1,rotation_error1]=cal_pose_error_noise(res,R,t);
            method_list(k).r1(:,i)= rotation_error1; 
            method_list(k).t1(:,i)= position_error1;
        end               
        showpercent(i,iterations);
    end
    %Now compute the mean of the error for each algorithm;
    for k=1:length(method_list)
        method_list(k).mean_r1(:,n)=mean(method_list(k).r1,2);
        method_list(k).mean_t1(:,n)=mean(method_list(k).t1,2);
    end
    fprintf('\n');
end


figure(1);
hold on;box on;
for k=1:length(method_list)
    plot(noise_levels,method_list(k).mean_r1(1,:),'color',method_list(k).color,'LineWidth',2);    
end
set(gca,'FontSize',14);
xlabel('Noise Level (Pixel)');ylabel('Rotation Error (Deg)');
legend(name,'Location','NorthWest');

figure(2);
hold on;box on;
for k=1:length(method_list)
    plot(noise_levels,method_list(k).mean_t1(1,:),'color',method_list(k).color,'LineWidth',2);    
end
set(gca,'FontSize',14);
xlabel('Noise Level (Pixel)');ylabel('Translation Error (%)');
legend(name,'Location','NorthWest');


