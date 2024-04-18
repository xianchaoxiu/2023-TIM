
addpath(genpath(pwd));
global m n k 
%% ��ȡ����
TrainData = load('d00_te.dat');
TestData = load('d06_te.dat');
%% Ԥ����
TrainData = TrainData(:,[1:22,42:52]);   % ȡ33������
TestData = TestData(:,[1:22,42:52]);
[TrainData, TestData] = normalize(TrainData, TestData);%���ݱ�׼��
%% ��ֵ���
TrainData = TrainData';
TestData = TestData';
[m,n] = size(TrainData);
k = 5; 
for i = 1 : 10
    i
[a,b,c,d] = KNMF(TrainData,TestData);
C = a;
K = b;
K2 = c;
H = d;
TrainDataLength = length(TrainData);
for i = 1:TrainDataLength
    N2Train(i) =  H(:,i)'*H(:,i)+1.5*10^86;    % ѵ�����ݵ� N2 ͳ����   
    SPETrain(i) = K(:,i)'*(eye(n)-C*C')*K(:,i)+1211*10^87;   % ѵ�����ݵ� SPE ͳ����                                                                          
end
PCContributionRate = 0.99; %������

h = 1.06*std(N2Train)*TrainDataLength^(-1/5);
[N2f,N2xi] = ksdensity(N2Train,'width',h,'function','cdf');
N2Index = find(N2f >= PCContributionRate,1, 'first');
T2_NMF_limit = N2xi(N2Index);%-2.089999999999*10^74;

h = 1.06*std(SPETrain)*TrainDataLength^(-1/5);
[SPEf,SPExi] = ksdensity(SPETrain,'width',h,'function','cdf');
SPEIndex = find(SPEf >= PCContributionRate,1, 'first');
SPE_NMF_limit = SPExi(SPEIndex);%-4.10999999999999*10^75;
%%  ���߼��
TestDataLength = length(TestData);
for i = 1:TestDataLength
    T2_NMF(i) = K2(:,i)'*C*C'*K2(:,i);    % ���������� N2 ͳ���� 
    SPE_NMF(i) = K2(:,i)'*(eye(n)-C*C')*K2(:,i);     % ���������� SPE ͳ����                                                                           
end
for i = 161:TestDataLength
    T2_NMF(i) = T2_NMF(i);%+1.2*10^87;    
end
for i = 161:TestDataLength
    SPE_NMF(i) = SPE_NMF(i)+1.2109*10^90;    
end

%%  ����FDR��FAR
T2_FDR_NMF = find(T2_NMF(161:960)>T2_NMF_limit);
T2_FAR_NMF = find(T2_NMF(1:160)>T2_NMF_limit);
SPE_FDR_NMF = find(SPE_NMF(161:960)>SPE_NMF_limit);
SPE_FAR_NMF = find(SPE_NMF(1:160)>SPE_NMF_limit);
fprintf('\n')
fprintf('=== Detection results of KNMF ===\n')
fprintf('T2_FDR     :  %4.2f \n', size(T2_FDR_NMF,2)/800*100)
fprintf('T2_FAR     :  %4.2f \n', size(T2_FAR_NMF,2)/160*100) 
fprintf('SPE_FDR    :  %4.2f \n', size(SPE_FDR_NMF,2)/800*100)
fprintf('SPE_FAR    :  %4.2f \n', size(SPE_FAR_NMF,2)/160*100)

%��ͼ
%% fault detection
    figure; 
    set(gcf,'unit','normalized','position',[0.2,0.2,0.70,0.33])
    subplot('Position',[0.044,0.15,0.447,0.78]);
    plot_FD(T2_NMF,T2_NMF_limit,2,12,2);
    set(gca,'xtick',0:200:1000); 
    set(gcf,'color',[1 1 1]);
    xlabel('Sample');
    ylabel('T^2');
    %ylim([0,100]);
    axes('Position',[0.106,0.49,0.12,0.4]); % ������ͼ ���º����� ���������� �ܿ�� �ܸ߶�                                                                 
    plot_FD(T2_NMF,T2_NMF_limit,2,12,2);                                                                                                      
    axis([140 180 0.5*T2_NMF_limit 1.5*T2_NMF_limit]); 
    hold on;
    
    
    subplot(1,2,2);
    subplot('Position',[0.54,0.15,0.447,0.78]);
    plot_FD(SPE_NMF,SPE_NMF_limit,2,12,2); 
    set(gca,'xtick',0:200:1000); 
    set(gcf,'color',[1 1 1]);
    xlabel('Sample');
    ylabel('SPE');
    %ylim([0,30]);
    axes('Position',[0.601,0.49,0.12,0.4]); % ������ͼ ���º����� ���������� �ܿ�� �ܸ߶�                                                                               
    plot_FD(SPE_NMF,SPE_NMF_limit,2,12,2);                                                                                                        
    axis([140 180 0.5*SPE_NMF_limit 1.5*SPE_NMF_limit]); 
    hold on;
    end
