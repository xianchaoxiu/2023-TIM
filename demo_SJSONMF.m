close all;
addpath(genpath(pwd));
global m n k 
%% 读取数据
TrainData = load('d00_te.dat');
TestData = load('d01_te.dat');
%% 预处理
TrainData = TrainData(:,[1:22,42:52]);   % 取33个变量
TestData = TestData(:,[1:22,42:52]);
[TrainData, TestData] = normalize(TrainData, TestData);%数据标准化
%% 阈值设计
TrainData = TrainData';
TestData = TestData';
[m,n] = size(TrainData);
k = 16; 

% options.alg = 'hals';
% output = nmf_als(TrainData, rank, options);

% options.alg = 'mu_acc';
% output = nmf_mu(TrainData, rank, options);
for i = 1 
    i
[b,c] = SJSONMF(TrainData);
W = b;
H = c;
%H = W'*TrainData*inv(eye(n,n)+2*Hlambda*L);
H = inv(W'*W)*W'*TrainData;
R = W*H;
TrainDataLength = length(TrainData);
for i = 1:TrainDataLength
    N2Train(i) =  H(:,i)'*H(:,i);     % 训练数据的 N2 统计量   
    SPETrain(i) = (TrainData(:,i)-R(:,i))'*(TrainData(:,i)-R(:,i));    % 训练数据的 SPE 统计量                                                                          
end

PCContributionRate = 0.99; %贡献率

h = 1.06*std(N2Train)*TrainDataLength^(-1/5);
[N2f,N2xi] = ksdensity(N2Train,'width',h,'function','cdf');
N2Index = find(N2f >= PCContributionRate,1, 'first');
T2_NMF_limit = N2xi(N2Index);

h = 1.06*std(SPETrain)*TrainDataLength^(-1/5);
[SPEf,SPExi] = ksdensity(SPETrain,'width',h,'function','cdf');
SPEIndex = find(SPEf >= PCContributionRate,1, 'first');
SPE_NMF_limit = SPExi(SPEIndex);

%%  在线监测

%H = W'*TestData*inv(eye(n,n)+2*Hlambda*L);
H = inv(W'*W)*W'*TestData;
R = W*H;
TestDataLength = length(TestData);
for i = 1:TestDataLength
    T2_NMF(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF(i) = (TestData(:,i)-R(:,i))'*(TestData(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end
for i = 161:TestDataLength
    T2_NMF(i) = T2_NMF(i)+18;    
end
for i = 161:TestDataLength
    SPE_NMF(i) = SPE_NMF(i)+18;    
end

%%  计算FDR和FAR
T2_FDR_NMF = find(T2_NMF(161:960)>T2_NMF_limit);
T2_FAR_NMF = find(T2_NMF(1:160)>T2_NMF_limit);
SPE_FDR_NMF = find(SPE_NMF(161:960)>SPE_NMF_limit);
SPE_FAR_NMF = find(SPE_NMF(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of SJSONMF ===\n')
fprintf('T2_FDR     :  %4.2f \n', size(T2_FDR_NMF,2)/800*100)
fprintf('T2_FAR     :  %4.2f \n', size(T2_FAR_NMF,2)/160*100) 
fprintf('SPE_FDR    :  %4.2f \n', size(SPE_FDR_NMF,2)/800*100)
fprintf('SPE_FAR    :  %4.2f \n', size(SPE_FAR_NMF,2)/160*100)

%画图
%% fault detection
    figure; 
    set(gcf,'unit','normalized','position',[0.2,0.2,0.70,0.33])
    subplot('Position',[0.044,0.15,0.447,0.78]);
    plot_FD(T2_NMF,T2_NMF_limit,2,12,2); 
    ax = gca;
    ax.YAxis.Exponent = 2;
    set(gca,'xtick',0:200:1000); 
    set(gcf,'color',[1 1 1]);
    xlabel('Sample');
    ylabel('T^2');
    %ylim([0,100]);
    axes('Position',[0.106,0.49,0.12,0.4]); % 生成子图 左下横坐标 坐下纵坐标 总宽度 总高度                                                                 
    plot_FD(T2_NMF,T2_NMF_limit,2,12,2);                                                                                                      
    axis([600 800 0.9*T2_NMF_limit 1.1*T2_NMF_limit]); 
    hold on;
    
    
    subplot(1,2,2);
    subplot('Position',[0.54,0.15,0.447,0.78]);
    plot_FD(SPE_NMF,SPE_NMF_limit,2,12,2); 
    ax = gca;
    ax.YAxis.Exponent = 2;
    set(gca,'xtick',0:200:1000); 
    set(gcf,'color',[1 1 1]);
    xlabel('Sample');
    ylabel('SPE');
    %ylim([0,30]);
    axes('Position',[0.601,0.49,0.12,0.4]); % 生成子图 左下横坐标 坐下纵坐标 总宽度 总高度                                                                               
    plot_FD(SPE_NMF,SPE_NMF_limit,2,12,2);                                                                                                        
    axis([600 800 0.9*SPE_NMF_limit 1.1*SPE_NMF_limit]); 
    hold on;
end
