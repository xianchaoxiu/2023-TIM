clc;
clear all; 
close all;
addpath(genpath(pwd));

%% 读取数据
TrainData = load('d00_te.dat'); 
TestData = load('d06_te.dat');

%% 预处理
TrainData = TrainData(:,[1:22,42:52]);   % 取33个变量
TestData = TestData(:,[1:22,42:52]);
[TrainData, TestData] = normalize(TrainData, TestData);  %数据标准化;
[X_row,X_col] = size(TrainData); 

%% 阈值设计
tic
[U,S,V] = svd(TrainData./sqrt(size(TrainData,1)-1));  
toc
D = S.^2;
size(V);
lambda = diag(D);
num_pc = 1; 
while sum(lambda(1:num_pc))/sum(lambda) < 0.85
    num_pc = num_pc+1;
end
D = diag(lambda(1:num_pc));
T2_PCA_limit=num_pc*(X_row-1)*(X_row+1)*finv(0.99,num_pc,X_row - num_pc)/(X_row*(X_row - num_pc)); 

for i = 1:3       
    theta(i) = sum((lambda(num_pc+1:X_col)).^i);
end
h0 = 1 - 2*theta(1)*theta(3)/(3*theta(2)^2);
ca = norminv(0.99,0,1); 
SPE_PCA_limit = theta(1)*(h0*ca*sqrt(2*theta(2))/theta(1) + 1 + theta(2)*h0*(h0 - 1)/theta(1)^2)^(1/h0); 

%% 在线监测
lambda = diag(lambda);
D = lambda(1:num_pc,1:num_pc);

P1 = V(:,1:num_pc);
[r,y] = size(P1*P1');
I = eye(r,y);

T2_PCA = zeros(X_row,1);
SPE_PCA = zeros(X_row,1);
for i = 1:X_row
    T2_PCA(i) = TestData(i,:)*P1*inv(D)*P1'*TestData(i,:)';                                      
    SPE_PCA(i) = TestData(i,:)*(I - P1*P1')*TestData(i,:)';                                                                                   
end

%%  计算FDR和FAR
T2_FDR_PCA = find(T2_PCA(161:960)>T2_PCA_limit);
T2_FAR_PCA = find(T2_PCA(1:160)>T2_PCA_limit);
SPE_FDR_PCA = find(SPE_PCA(161:960)>SPE_PCA_limit);
SPE_FAR_PCA = find(SPE_PCA(1:160)>SPE_PCA_limit);
fprintf('\n')
fprintf('=== Detection results of PCA ===\n')
fprintf('T2_FDR     :  %4.2f \n', size(T2_FDR_PCA,1)/800*100)
fprintf('T2_FAR     :  %4.2f \n', size(T2_FAR_PCA,1)/160*100) 
fprintf('SPE_FDR    :  %4.2f \n', size(SPE_FDR_PCA,1)/800*100)
fprintf('SPE_FAR    :  %4.2f \n', size(SPE_FAR_PCA,1)/160*100)

%% fault detection
    figure; 
    set(gcf,'unit','normalized','position',[0.2,0.2,0.70,0.33])
    subplot('Position',[0.044,0.15,0.447,0.78]);
    plot_FD(T2_PCA,T2_PCA_limit,2,12,2); 
    ax = gca;
    ax.YAxis.Exponent = 3;
    set(gca,'xtick',0:200:1000); 
    set(gcf,'color',[1 1 1]);
    xlabel('Sample');
    ylabel('T^2');
    %ylim([0,100]);
    axes('Position',[0.106,0.49,0.12,0.4]); % 生成子图 左下横坐标 坐下纵坐标 总宽度 总高度                                                                 
    plot_FD(T2_PCA,T2_PCA_limit,2,12,2);                                                                                                      
    axis([140 180 T2_PCA_limit-10 T2_PCA_limit+10]); 
    hold on;
    
    
    subplot(1,2,2);
    subplot('Position',[0.54,0.15,0.447,0.78]);
    plot_FD(SPE_PCA,SPE_PCA_limit,2,12,2); 
    set(gca,'FontSize',12.6)
    set(gca,'xtick',0:200:1000); 
    set(gcf,'color',[1 1 1]);
    xlabel('Sample');
    ylabel('SPE');
    %ylim([0,30]);
    axes('Position',[0.601,0.49,0.12,0.4]); % 生成子图 左下横坐标 坐下纵坐标 总宽度 总高度                                                                               
    plot_FD(SPE_PCA,SPE_PCA_limit,2,12,2);                                                                                                        
    axis([140 180 SPE_PCA_limit-10 SPE_PCA_limit+10]); 
    hold on;

    

  