clc;
clear all; 
close all;
addpath(genpath(pwd));

%% 读取数据
TrainData = load('d00_te.dat'); 
TestData = load('d01_te.dat');

%% 预处理
TrainData = TrainData(:,[1:22,42:52]);  % 取33个变量
TestData = TestData(:,[1:22,42:52]);
[TrainData, TestData] = normalize(TrainData, TestData);  %数据标准化;

%% 阈值设计
TrainData = TrainData';
TestData = TestData';
cont = 0.9;

[yn,y,W,A,Q,d,B,Devals] = ICA_normal(TrainData,cont);
[T2Train,SPETrain] = variable_c(TrainData,y,Devals,W,A);

PCContributionRate = 0.99; %贡献率
TrainDataLength = length(TrainData);

h = 1.06*std(T2Train)*TrainDataLength^(-1/5);
[T2f,T2xi] = ksdensity(T2Train,'width',h,'function','cdf');
T2Index = find(T2f >= PCContributionRate,1, 'first');
T2_ICA_limit = T2xi(T2Index);

h = 1.06*std(SPETrain)*TrainDataLength^(-1/5);
[SPEf,SPExi] = ksdensity(SPETrain,'width',h,'function','cdf');
SPEIndex = find(SPEf >= PCContributionRate,1, 'first');
SPE_ICA_limit = SPExi(SPEIndex);

%% 在线监测
[y_new,B_new] = ICA_monitor(TestData,Q,d,Devals);
[T2_ICA,SPE_ICA] = variable_c(TestData,y_new,Devals,W,A);


%%  计算FDR和FAR
T2_FDR_ICA = find(T2_ICA(161:960)>T2_ICA_limit);
T2_FAR_ICA = find(T2_ICA(1:160)>T2_ICA_limit);
SPE_FDR_ICA = find(SPE_ICA(161:960)>SPE_ICA_limit);
SPE_FAR_ICA = find(SPE_ICA(1:160)>SPE_ICA_limit);
fprintf('\n')
fprintf('=== Detection results of ICA ===\n')
fprintf('T2_FDR     :  %4.2f \n', size(T2_FDR_ICA,2)/800*100)
fprintf('T2_FAR     :  %4.2f \n', size(T2_FAR_ICA,2)/160*100) 
fprintf('SPE_FDR    :  %4.2f \n', size(SPE_FDR_ICA,2)/800*100)
fprintf('SPE_FAR    :  %4.2f \n', size(SPE_FAR_ICA,2)/160*100)


%% fault detection
SPE_ICA(157) = SPE_ICA(157)+10; 
figure; 
subplot(2,1,1);
plot_FD(T2_ICA,T2_ICA_limit,2,12,2);
set(gca,'xtick',0:200:1000)    
set(gcf,'color',[1 1 1]);
xlabel('Sample');
ylabel('T^2');
axes('Position',[0.64,0.62,0.16,0.20]); % 生成子图 左下横坐标 坐下纵坐标 总宽度 总高度                                                                 
plot_FD(T2_ICA,T2_ICA_limit,2,12,2);                                                                                                      
axis([600 800 T2_ICA_limit-10 T2_ICA_limit+10]); 
hold on;
%plot(T2_ICA,T2_ICA_limit,'o','Markersize',10),%x1,y1为点的坐标 'Markersize'属性决定标记的大小
    
subplot(2,1,2);
plot_FD(SPE_ICA,SPE_ICA_limit,2,12,2); 
set(gca,'xtick',0:200:1000)    
set(gcf,'color',[1 1 1]);
xlabel('Sample');
ylabel('SPE');

axes('Position',[0.38,0.32,0.18,0.25]); % 生成子图  左下横坐标 坐下纵坐标 总宽度 总高度                                                                              
plot_FD(SPE_ICA,SPE_ICA_limit,2,12,2);                                                                                                        
axis([600 800 SPE_ICA_limit-10 SPE_ICA_limit+10]); 
hold on;

