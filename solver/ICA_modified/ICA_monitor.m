function [y,B]=ICA_monitor(X,Q,d,Devals)
[m,n]=size(X); %获取输入矩阵的行/列数，行数为观测数据的数目，列数为采样点数

z=Q*X;   %正交矩阵

%---------正交矩阵B------
B=[eye(d,d);zeros(m-d,d)];% 初始化列向量B的寄存矩阵,B=[b1  b2  ...   bd]
for r=1:d                            % 迭代求取每一个独立元
    b=B(:,r);
    i=1;maxIterationsNum=500;j=1;   % 设置最大迭代次数（即对于每个独立分量而言迭代均不超过此次数）
    while i<=maxIterationsNum+1
        if i==maxIterationsNum     % 循环结束处理
            fprintf('\n第%d分量在%d次迭代内并不收敛。',r,maxIterationsNum);
            break;
        end
        b0=b;                     
        t=z'*b;
        g=tanh(t);
        dg=1-tanh(t).^2;
        b=z*g/n-mean(dg)*b; % 核心公式，见参考文献p5
        %%%%%%%%%%%% 对b正交化
        for k=1:r-1
            b=b-(b'*B(:,k))*B(:,k);
        end                    
        b=b/norm(b);
        %%%%%%%%%%%  收敛性判断
        if norm(b0-b)/norm(b0)<10^-5 
           % abs(abs(b'*b0)-1)<10^-5        % 如果收敛，则保存b
               B(:,r)=b; 
             break;
        end
        i=i+1;        
    end 
end

yn=B'*z; % ICA分离的独立元
y=sqrt(Devals)*yn;