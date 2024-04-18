clc;
clear all; 
close all;
addpath(genpath(pwd));
global m n k 
%% 读取数据
TrainData = load('d00_te.dat'); 
TestData1 = load('d01_te.dat');
TestData2 = load('d02_te.dat');
TestData3 = load('d03_te.dat');
TestData4 = load('d04_te.dat');
TestData5 = load('d05_te.dat');
TestData6 = load('d06_te.dat');
TestData7 = load('d07_te.dat');
TestData8 = load('d08_te.dat');
TestData9 = load('d09_te.dat');
TestData10 = load('d10_te.dat');
TestData11 = load('d11_te.dat');
TestData12 = load('d12_te.dat');
TestData13 = load('d13_te.dat');
TestData14 = load('d14_te.dat');
TestData15 = load('d15_te.dat');
TestData16 = load('d16_te.dat');
TestData17 = load('d17_te.dat');
TestData18 = load('d18_te.dat');
TestData19 = load('d19_te.dat');
TestData20 = load('d20_te.dat');
TestData21 = load('d21_te.dat');

%% 预处理
TrainData = TrainData(:,[1:22,42:52]);   % 取33个变量
TestData1 = TestData1(:,[1:22,42:52]);
TestData2 = TestData2(:,[1:22,42:52]);
TestData3 = TestData3(:,[1:22,42:52]);
TestData4 = TestData4(:,[1:22,42:52]);
TestData5 = TestData5(:,[1:22,42:52]);
TestData6 = TestData6(:,[1:22,42:52]);
TestData7 = TestData7(:,[1:22,42:52]);
TestData8 = TestData8(:,[1:22,42:52]);
TestData9 = TestData9(:,[1:22,42:52]);
TestData10 = TestData10(:,[1:22,42:52]);
TestData11 = TestData11(:,[1:22,42:52]);
TestData12 = TestData12(:,[1:22,42:52]);
TestData13 = TestData13(:,[1:22,42:52]);
TestData14 = TestData14(:,[1:22,42:52]);
TestData15 = TestData15(:,[1:22,42:52]);
TestData16 = TestData16(:,[1:22,42:52]);
TestData17 = TestData17(:,[1:22,42:52]);
TestData18 = TestData18(:,[1:22,42:52]);
TestData19 = TestData19(:,[1:22,42:52]);
TestData20 = TestData20(:,[1:22,42:52]);
TestData21 = TestData21(:,[1:22,42:52]);
% [TrainData, TestData] = normalize(TrainData, TestData);  %数据标准化

%% 阈值设计
TrainData = TrainData';
TestData1 = TestData1';
TestData2 = TestData2';
TestData3 = TestData3';
TestData4 = TestData4';
TestData5 = TestData5';
TestData6 = TestData6';
TestData7 = TestData7';
TestData8 = TestData8';
TestData9 = TestData9';
TestData10 = TestData10';
TestData11 = TestData11';
TestData12 = TestData12';
TestData13 = TestData13';
TestData14 = TestData14';
TestData15 = TestData15';
TestData16 = TestData16';
TestData17 = TestData17';
TestData18 = TestData18';
TestData19 = TestData19';
TestData20 = TestData20';
TestData21 = TestData21';
rank = 17;

% options.alg = 'hals';
% output = nmf_als(TrainData, rank, options);

% options.alg = 'mu_acc';
% output = nmf_mu(TrainData, rank, options);
for i=1 
    tic
options.alg = 'nmf_hals_so';
[x, infos] = nmf_hals_so(TrainData, rank, options);
toc
%[U_final, V_final, nIter_final, objhistory_final] = GNMF(TrainData, rank,options);
W = x.W;
H = pinv(W'*W)*W'*TrainData;
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
H = pinv(W'*W)*W'*TestData1;
R = W*H;
TestDataLength = length(TestData1);
for i = 1:TestDataLength
    T2_NMF1(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF1(i) = (TestData1(:,i)-R(:,i))'*(TestData1(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF1 = find(T2_NMF1(161:960)>=T2_NMF_limit);
T2_FAR_NMF1 = find(T2_NMF1(1:160)>T2_NMF_limit);
SPE_FDR_NMF1 = find(SPE_NMF1(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF1 = find(SPE_NMF1(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR1     :  %4.2f \n', size(T2_FDR_NMF1,2)/800*100)
fprintf('T2_FAR1     :  %4.2f \n', size(T2_FAR_NMF1,2)/160*100) 
fprintf('SPE_FDR1    :  %4.2f \n', size(SPE_FDR_NMF1,2)/800*100)
fprintf('SPE_FAR1    :  %4.2f \n', size(SPE_FAR_NMF1,2)/160*100)

%%  在线监测
H = pinv(W'*W)*W'*TestData2;
R = W*H;
TestDataLength = length(TestData2);
for i = 1:TestDataLength
    T2_NMF2(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF2(i) = (TestData2(:,i)-R(:,i))'*(TestData2(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF2 = find(T2_NMF2(161:960)>=T2_NMF_limit);
T2_FAR_NMF2 = find(T2_NMF2(1:160)>T2_NMF_limit);
SPE_FDR_NMF2 = find(SPE_NMF2(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF2 = find(SPE_NMF2(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR2     :  %4.2f \n', size(T2_FDR_NMF2,2)/800*100)
fprintf('T2_FAR2     :  %4.2f \n', size(T2_FAR_NMF2,2)/160*100) 
fprintf('SPE_FDR2    :  %4.2f \n', size(SPE_FDR_NMF2,2)/800*100)
fprintf('SPE_FAR2    :  %4.2f \n', size(SPE_FAR_NMF2,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData3;
R = W*H;
TestDataLength = length(TestData3);
for i = 1:TestDataLength
    T2_NMF3(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF3(i) = (TestData3(:,i)-R(:,i))'*(TestData3(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF3 = find(T2_NMF3(161:960)>=T2_NMF_limit);
T2_FAR_NMF3 = find(T2_NMF3(1:160)>T2_NMF_limit);
SPE_FDR_NMF3 = find(SPE_NMF3(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF3 = find(SPE_NMF3(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR3     :  %4.2f \n', size(T2_FDR_NMF3,2)/800*100)
fprintf('T2_FAR3     :  %4.2f \n', size(T2_FAR_NMF3,2)/160*100) 
fprintf('SPE_FDR3    :  %4.2f \n', size(SPE_FDR_NMF3,2)/800*100)
fprintf('SPE_FAR3    :  %4.2f \n', size(SPE_FAR_NMF3,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData4;
R = W*H;
TestDataLength = length(TestData4);
for i = 1:TestDataLength
    T2_NMF4(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF4(i) = (TestData4(:,i)-R(:,i))'*(TestData4(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF4 = find(T2_NMF4(161:960)>=T2_NMF_limit);
T2_FAR_NMF4 = find(T2_NMF4(1:160)>T2_NMF_limit);
SPE_FDR_NMF4 = find(SPE_NMF4(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF4 = find(SPE_NMF4(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR4     :  %4.2f \n', size(T2_FDR_NMF4,2)/800*100)
fprintf('T2_FAR4    :  %4.2f \n', size(T2_FAR_NMF4,2)/160*100) 
fprintf('SPE_FDR4    :  %4.2f \n', size(SPE_FDR_NMF4,2)/800*100)
fprintf('SPE_FAR4    :  %4.2f \n', size(SPE_FAR_NMF4,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData5;
R = W*H;
TestDataLength = length(TestData5);
for i = 1:TestDataLength
    T2_NMF5(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF5(i) = (TestData5(:,i)-R(:,i))'*(TestData5(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF5 = find(T2_NMF5(161:960)>=T2_NMF_limit);
T2_FAR_NMF5 = find(T2_NMF5(1:160)>T2_NMF_limit);
SPE_FDR_NMF5 = find(SPE_NMF5(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF5 = find(SPE_NMF5(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR5     :  %4.2f \n', size(T2_FDR_NMF5,2)/800*100)
fprintf('T2_FAR5     :  %4.2f \n', size(T2_FAR_NMF5,2)/160*100) 
fprintf('SPE_FDR5    :  %4.2f \n', size(SPE_FDR_NMF5,2)/800*100)
fprintf('SPE_FAR5    :  %4.2f \n', size(SPE_FAR_NMF5,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData6;
R = W*H;
TestDataLength = length(TestData6);
for i = 1:TestDataLength
    T2_NMF6(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF6(i) = (TestData6(:,i)-R(:,i))'*(TestData6(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF6 = find(T2_NMF6(161:960)>=T2_NMF_limit);
T2_FAR_NMF6 = find(T2_NMF6(1:160)>T2_NMF_limit);
SPE_FDR_NMF6 = find(SPE_NMF6(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF6 = find(SPE_NMF6(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR6     :  %4.2f \n', size(T2_FDR_NMF6,2)/800*100)
fprintf('T2_FAR6     :  %4.2f \n', size(T2_FAR_NMF6,2)/160*100) 
fprintf('SPE_FDR6    :  %4.2f \n', size(SPE_FDR_NMF6,2)/800*100)
fprintf('SPE_FAR6    :  %4.2f \n', size(SPE_FAR_NMF6,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData7;
R = W*H;
TestDataLength = length(TestData7);
for i = 1:TestDataLength
    T2_NMF7(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF7(i) = (TestData7(:,i)-R(:,i))'*(TestData7(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF7 = find(T2_NMF7(161:960)>=T2_NMF_limit);
T2_FAR_NMF7 = find(T2_NMF7(1:160)>T2_NMF_limit);
SPE_FDR_NMF7 = find(SPE_NMF7(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF7 = find(SPE_NMF7(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR7     :  %4.2f \n', size(T2_FDR_NMF7,2)/800*100)
fprintf('T2_FAR7     :  %4.2f \n', size(T2_FAR_NMF7,2)/160*100) 
fprintf('SPE_FDR7    :  %4.2f \n', size(SPE_FDR_NMF7,2)/800*100)
fprintf('SPE_FAR7    :  %4.2f \n', size(SPE_FAR_NMF7,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData8;
R = W*H;
TestDataLength = length(TestData8);
for i = 1:TestDataLength
    T2_NMF8(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF8(i) = (TestData8(:,i)-R(:,i))'*(TestData8(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF8 = find(T2_NMF8(161:960)>=T2_NMF_limit);
T2_FAR_NMF8 = find(T2_NMF8(1:160)>T2_NMF_limit);
SPE_FDR_NMF8 = find(SPE_NMF8(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF8 = find(SPE_NMF8(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR8     :  %4.2f \n', size(T2_FDR_NMF8,2)/800*100)
fprintf('T2_FAR8     :  %4.2f \n', size(T2_FAR_NMF8,2)/160*100) 
fprintf('SPE_FDR8    :  %4.2f \n', size(SPE_FDR_NMF8,2)/800*100)
fprintf('SPE_FAR8    :  %4.2f \n', size(SPE_FAR_NMF8,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData9;
R = W*H;
TestDataLength = length(TestData9);
for i = 1:TestDataLength
    T2_NMF9(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF9(i) = (TestData9(:,i)-R(:,i))'*(TestData9(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF9 = find(T2_NMF9(161:960)>=T2_NMF_limit);
T2_FAR_NMF9 = find(T2_NMF9(1:160)>T2_NMF_limit);
SPE_FDR_NMF9 = find(SPE_NMF9(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF9 = find(SPE_NMF9(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR9     :  %4.2f \n', size(T2_FDR_NMF9,2)/800*100)
fprintf('T2_FAR9     :  %4.2f \n', size(T2_FAR_NMF9,2)/160*100) 
fprintf('SPE_FDR9    :  %4.2f \n', size(SPE_FDR_NMF9,2)/800*100)
fprintf('SPE_FAR9    :  %4.2f \n', size(SPE_FAR_NMF9,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData10;
R = W*H;
TestDataLength = length(TestData10);
for i = 1:TestDataLength
    T2_NMF10(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF10(i) = (TestData10(:,i)-R(:,i))'*(TestData10(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF10 = find(T2_NMF10(161:960)>=T2_NMF_limit);
T2_FAR_NMF10 = find(T2_NMF10(1:160)>T2_NMF_limit);
SPE_FDR_NMF10 = find(SPE_NMF10(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF10 = find(SPE_NMF10(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR10     :  %4.2f \n', size(T2_FDR_NMF10,2)/800*100)
fprintf('T2_FAR10     :  %4.2f \n', size(T2_FAR_NMF10,2)/160*100) 
fprintf('SPE_FDR10    :  %4.2f \n', size(SPE_FDR_NMF10,2)/800*100)
fprintf('SPE_FAR10    :  %4.2f \n', size(SPE_FAR_NMF10,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData11;
R = W*H;
TestDataLength = length(TestData11);
for i = 1:TestDataLength
    T2_NMF11(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF11(i) = (TestData11(:,i)-R(:,i))'*(TestData11(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF11 = find(T2_NMF11(161:960)>=T2_NMF_limit);
T2_FAR_NMF11 = find(T2_NMF11(1:160)>T2_NMF_limit);
SPE_FDR_NMF11 = find(SPE_NMF11(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF11 = find(SPE_NMF11(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR11     :  %4.2f \n', size(T2_FDR_NMF11,2)/800*100)
fprintf('T2_FAR11     :  %4.2f \n', size(T2_FAR_NMF11,2)/160*100) 
fprintf('SPE_FDR11    :  %4.2f \n', size(SPE_FDR_NMF11,2)/800*100)
fprintf('SPE_FAR11    :  %4.2f \n', size(SPE_FAR_NMF11,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData12;
R = W*H;
TestDataLength = length(TestData12);
for i = 1:TestDataLength
    T2_NMF12(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF12(i) = (TestData12(:,i)-R(:,i))'*(TestData12(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF12 = find(T2_NMF12(161:960)>=T2_NMF_limit);
T2_FAR_NMF12 = find(T2_NMF12(1:160)>T2_NMF_limit);
SPE_FDR_NMF12 = find(SPE_NMF12(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF12 = find(SPE_NMF12(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR12     :  %4.2f \n', size(T2_FDR_NMF12,2)/800*100)
fprintf('T2_FAR12     :  %4.2f \n', size(T2_FAR_NMF12,2)/160*100) 
fprintf('SPE_FDR12    :  %4.2f \n', size(SPE_FDR_NMF12,2)/800*100)
fprintf('SPE_FAR12    :  %4.2f \n', size(SPE_FAR_NMF12,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData13;
R = W*H;
TestDataLength = length(TestData13);
for i = 1:TestDataLength
    T2_NMF13(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF13(i) = (TestData13(:,i)-R(:,i))'*(TestData13(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF13 = find(T2_NMF13(161:960)>=T2_NMF_limit);
T2_FAR_NMF13 = find(T2_NMF13(1:160)>T2_NMF_limit);
SPE_FDR_NMF13 = find(SPE_NMF13(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF13 = find(SPE_NMF13(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR13     :  %4.2f \n', size(T2_FDR_NMF13,2)/800*100)
fprintf('T2_FAR13     :  %4.2f \n', size(T2_FAR_NMF13,2)/160*100) 
fprintf('SPE_FDR13    :  %4.2f \n', size(SPE_FDR_NMF13,2)/800*100)
fprintf('SPE_FAR13    :  %4.2f \n', size(SPE_FAR_NMF13,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData14;
R = W*H;
TestDataLength = length(TestData14);
for i = 1:TestDataLength
    T2_NMF14(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF14(i) = (TestData14(:,i)-R(:,i))'*(TestData14(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF14 = find(T2_NMF14(161:960)>=T2_NMF_limit);
T2_FAR_NMF14 = find(T2_NMF14(1:160)>T2_NMF_limit);
SPE_FDR_NMF14 = find(SPE_NMF14(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF14 = find(SPE_NMF14(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR14     :  %4.2f \n', size(T2_FDR_NMF14,2)/800*100)
fprintf('T2_FAR14    :  %4.2f \n', size(T2_FAR_NMF14,2)/160*100) 
fprintf('SPE_FDR14    :  %4.2f \n', size(SPE_FDR_NMF14,2)/800*100)
fprintf('SPE_FAR14    :  %4.2f \n', size(SPE_FAR_NMF14,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData15;
R = W*H;
TestDataLength = length(TestData15);
for i = 1:TestDataLength
    T2_NMF15(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF15(i) = (TestData15(:,i)-R(:,i))'*(TestData15(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF15 = find(T2_NMF15(161:960)>=T2_NMF_limit);
T2_FAR_NMF15 = find(T2_NMF15(1:160)>T2_NMF_limit);
SPE_FDR_NMF15 = find(SPE_NMF15(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF15 = find(SPE_NMF15(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR15     :  %4.2f \n', size(T2_FDR_NMF15,2)/800*100)
fprintf('T2_FAR15     :  %4.2f \n', size(T2_FAR_NMF15,2)/160*100) 
fprintf('SPE_FDR15    :  %4.2f \n', size(SPE_FDR_NMF15,2)/800*100)
fprintf('SPE_FAR15    :  %4.2f \n', size(SPE_FAR_NMF15,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData16;
R = W*H;
TestDataLength = length(TestData16);
for i = 1:TestDataLength
    T2_NMF16(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF16(i) = (TestData16(:,i)-R(:,i))'*(TestData16(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF16 = find(T2_NMF16(161:960)>=T2_NMF_limit);
T2_FAR_NMF16 = find(T2_NMF16(1:160)>T2_NMF_limit);
SPE_FDR_NMF16 = find(SPE_NMF16(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF16 = find(SPE_NMF16(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR16     :  %4.2f \n', size(T2_FDR_NMF16,2)/800*100)
fprintf('T2_FAR16     :  %4.2f \n', size(T2_FAR_NMF16,2)/160*100) 
fprintf('SPE_FDR16    :  %4.2f \n', size(SPE_FDR_NMF16,2)/800*100)
fprintf('SPE_FAR16    :  %4.2f \n', size(SPE_FAR_NMF16,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData17;
R = W*H;
TestDataLength = length(TestData17);
for i = 1:TestDataLength
    T2_NMF17(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF17(i) = (TestData17(:,i)-R(:,i))'*(TestData17(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF17 = find(T2_NMF17(161:960)>=T2_NMF_limit);
T2_FAR_NMF17 = find(T2_NMF17(1:160)>T2_NMF_limit);
SPE_FDR_NMF17 = find(SPE_NMF17(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF17 = find(SPE_NMF17(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR17     :  %4.2f \n', size(T2_FDR_NMF17,2)/800*100)
fprintf('T2_FAR17     :  %4.2f \n', size(T2_FAR_NMF17,2)/160*100) 
fprintf('SPE_FDR17    :  %4.2f \n', size(SPE_FDR_NMF17,2)/800*100)
fprintf('SPE_FAR17    :  %4.2f \n', size(SPE_FAR_NMF17,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData18;
R = W*H;
TestDataLength = length(TestData18);
for i = 1:TestDataLength
    T2_NMF18(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF18(i) = (TestData18(:,i)-R(:,i))'*(TestData18(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF18 = find(T2_NMF18(161:960)>=T2_NMF_limit);
T2_FAR_NMF18 = find(T2_NMF18(1:160)>T2_NMF_limit);
SPE_FDR_NMF18 = find(SPE_NMF18(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF18 = find(SPE_NMF18(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR18     :  %4.2f \n', size(T2_FDR_NMF18,2)/800*100)
fprintf('T2_FAR18     :  %4.2f \n', size(T2_FAR_NMF18,2)/160*100) 
fprintf('SPE_FDR18    :  %4.2f \n', size(SPE_FDR_NMF18,2)/800*100)
fprintf('SPE_FAR18    :  %4.2f \n', size(SPE_FAR_NMF18,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData19;
R = W*H;
TestDataLength = length(TestData19);
for i = 1:TestDataLength
    T2_NMF19(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF19(i) = (TestData19(:,i)-R(:,i))'*(TestData19(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF19 = find(T2_NMF19(161:960)>=T2_NMF_limit);
T2_FAR_NMF19 = find(T2_NMF19(1:160)>T2_NMF_limit);
SPE_FDR_NMF19 = find(SPE_NMF19(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF19 = find(SPE_NMF19(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR19     :  %4.2f \n', size(T2_FDR_NMF19,2)/800*100)
fprintf('T2_FAR19     :  %4.2f \n', size(T2_FAR_NMF19,2)/160*100) 
fprintf('SPE_FDR19    :  %4.2f \n', size(SPE_FDR_NMF19,2)/800*100)
fprintf('SPE_FAR19    :  %4.2f \n', size(SPE_FAR_NMF19,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData20;
R = W*H;
TestDataLength = length(TestData20);
for i = 1:TestDataLength
    T2_NMF20(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF20(i) = (TestData20(:,i)-R(:,i))'*(TestData20(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF20 = find(T2_NMF20(161:960)>=T2_NMF_limit);
T2_FAR_NMF20 = find(T2_NMF20(1:160)>T2_NMF_limit);
SPE_FDR_NMF20 = find(SPE_NMF20(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF20 = find(SPE_NMF20(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR20     :  %4.2f \n', size(T2_FDR_NMF20,2)/800*100)
fprintf('T2_FAR20     :  %4.2f \n', size(T2_FAR_NMF20,2)/160*100) 
fprintf('SPE_FDR20    :  %4.2f \n', size(SPE_FDR_NMF20,2)/800*100)
fprintf('SPE_FAR20    :  %4.2f \n', size(SPE_FAR_NMF20,2)/160*100)

%%  在线监测
H = inv(W'*W)*W'*TestData21;
R = W*H;
TestDataLength = length(TestData21);
for i = 1:TestDataLength
    T2_NMF21(i) =  H(:,i)'*H(:,i);     % 测试样本的 N2 统计量 
    SPE_NMF21(i) = (TestData21(:,i)-R(:,i))'*(TestData21(:,i)-R(:,i));     % 测试样本的 SPE 统计量                                                                           
end

%%  计算FDR和FAR
T2_FDR_NMF21 = find(T2_NMF21(161:960)>=T2_NMF_limit);
T2_FAR_NMF21 = find(T2_NMF21(1:160)>T2_NMF_limit);
SPE_FDR_NMF21 = find(SPE_NMF21(161:960)>=SPE_NMF_limit);
SPE_FAR_NMF21 = find(SPE_NMF21(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of GNMF ===\n')
fprintf('T2_FDR21     :  %4.2f \n', size(T2_FDR_NMF21,2)/800*100)
fprintf('T2_FAR21     :  %4.2f \n', size(T2_FAR_NMF21,2)/160*100) 
fprintf('SPE_FDR21    :  %4.2f \n', size(SPE_FDR_NMF21,2)/800*100)
fprintf('SPE_FAR21    :  %4.2f \n', size(SPE_FAR_NMF21,2)/160*100)


end
%   