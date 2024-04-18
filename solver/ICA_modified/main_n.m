clear all
clc
cont=0.9;
alpha=0.99;
c_aphla=2.32;    %对应于0.99，置信水平为95%对应分位点为1.645
%建模
xn=textread('d10.dat');
x=zscore(xn);
x=x';
[m,n]=size(x);
[yn,y,W,A,Q,d,B,Devals]=ICA_normal(x,cont);
[I2,SPE]=variable_c(x,y,Devals,W,A);
[f1,x1,u1]=ksdensity(I2);%I2单变量核密度估计
ConInt1=ComCon(f1,x1,alpha);
I2_limit=ConInt1(2);
SPE_limit=ksdensity(SPE,alpha,'function','icdf');
alpha=0.99;


%在线监控
Xn=textread('d10_te.dat');
X=zscore(Xn);
X=X';
fault_I2_num=[];
fault_SPE_num=[];
[y_new,B_new]=ICA_monitor(X,Q,d,Devals);
[I2_new,SPE_new]=variable_c(X,y_new,Devals,W,A);
for i=1:n
    if I2_new(i)>I2_limit
        fault_I2_num=[fault_I2_num,i];
    end;
    if SPE_new(i)>SPE_limit
        fault_SPE_num=[fault_SPE_num,i];
    end;
end;
fault_I2_num
fault_SPE_num

figure(1)
plot(1:length(I2_new),I2_new,'k*-');
hold on
plot(1:length(I2_new),ones(length(I2_new),1)*I2_limit,'r-');
xlabel('样本');ylabel('样本点的I2贡献值');
legend('I2统计量贡献','I2统计量控制限');
figure(2)
plot(1:length(SPE_new),SPE_new,'k*-');
hold on
plot(1:length(SPE_new),ones(length(SPE_new),1)*SPE_limit,'r-');
legend('SPE统计量贡献','SPE统计量控制限');
xlabel('样本');ylabel('样本点的SPE贡献值');

%故障诊断
new1=fault_I2_num(1);
new2=fault_SPE_num(1);
Xc_new=A*W*X;
cont_I2=yn(:,new1).^2;  
cont_SPE=(X(:,new2)-Xc_new(:,new2)).^2;
figure(3)
bar(cont_I2);
xlabel('变量');ylabel('变量对I2贡献');
figure(4)
bar(cont_SPE);
xlabel('变量');ylabel('变量对SPE贡献');