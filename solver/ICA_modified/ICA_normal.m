function [yn,y,W,A,Q,d,B,Devals]=ICA_normal(X,cont)
[m,n]=size(X); %��ȡ����������/����������Ϊ�۲����ݵ���Ŀ������Ϊ��������
%---------�׻�/��------
Cx = cov(X',1);    %����Э�������Cx 
[evecs,evals,V] = svd(Cx); %����Cx������ֵ����������

Q=evals^(-1/2)*evecs';   %�׻�����
z=Q*X;   %��������

Evals=diag(evals);
for i=1:m
    ICa(i)=sum(Evals(1:i))/sum(Evals);
    if ICa(i)>=cont
        d=i;
        break;
    end;
end;
Devals=diag(Evals(1:d));
Devecs=evecs(:,1:d);

%---------��������B------
B=[eye(d,d);zeros(m-d,d)];% ��ʼ��������B�ļĴ����,B=[b1  b2  ...   bd]
for r=1:d                            % ������ȡÿһ������Ԫ
    b=B(:,r);
    i=1;maxIterationsNum=500;j=1;   % ����������������������ÿ�������������Ե������������˴�����
    while i<=maxIterationsNum+1
        if i==maxIterationsNum     % ѭ����������
            fprintf('\n��%d������%d�ε����ڲ���������',r,maxIterationsNum);
            break;
        end
        b0=b;                     
        t=z'*b;
        g=tanh(t);
        dg=1-tanh(t).^2;
        b=z*g/n-mean(dg)*b; % ���Ĺ�ʽ�����ο�����p5
        %%%%%%%%%%%% ��b������
        for k=1:r-1
            b=b-(b'*B(:,k))*B(:,k);
        end                    
        b=b/norm(b);
        %%%%%%%%%%%  �������ж�
        if norm(b0-b)/norm(b0)<10^-5 
            % abs(abs(b'*b0)-1)<1e-6         % ����������򱣴�b
               B(:,r)=b; 
             break;
        end
        i=i+1;        
    end 
end

yn=B'*z; % ICA����Ķ���Ԫ
y=sqrt(Devals)*yn;
W=sqrt(Devals)*B'*Q;    % X=W*S
A=evecs*sqrt(evals)*B*Devals^(-1/2);