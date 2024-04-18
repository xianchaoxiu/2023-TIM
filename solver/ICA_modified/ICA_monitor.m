function [y,B]=ICA_monitor(X,Q,d,Devals)
[m,n]=size(X); %��ȡ����������/����������Ϊ�۲����ݵ���Ŀ������Ϊ��������

z=Q*X;   %��������

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
           % abs(abs(b'*b0)-1)<10^-5        % ����������򱣴�b
               B(:,r)=b; 
             break;
        end
        i=i+1;        
    end 
end

yn=B'*z; % ICA����Ķ���Ԫ
y=sqrt(Devals)*yn;