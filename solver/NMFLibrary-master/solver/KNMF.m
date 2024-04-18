function [C,K,K2,H] = KNMF(X,Y)
global n k 
d = 10;
r = 2000;
epslion = 10^-5;
C = rand(n,k);
H = rand(k,n);

% 核函数
K = zeros(n,n);
for i = 1 : n
    for j = 1 : n
        % K(i,j) = sqrt(1/(1+ norm(X(:,i)-X(:,j))^2));
        % K(i,j) = -log(norm(X(:,i)-X(:,j))^10+100); % 可调
        % K(i,j) = sqrt(norm(X(:,i)-X(:,j))^2+0.01); % 可调
        % K(i,j) = 1- norm(X(:,i)-X(:,j))^2/(norm(X(:,i)-X(:,j))^2+0.1);
        % K(i,j) = tanh(5*X(:,i)'*X(:,j)+1);
        K(i,j) = (r+50*X(:,i)'*X(:,j))^d; % 可调
    end
end

% 核函数
K2 = zeros(n,n);
for i = 1 : n
   for j = 1 : n
       % K2(i,j) = sqrt(1/(1+ norm(Y(:,i)-Y(:,j))^2));
       % K2(i,j) = -log(norm(Y(:,i)-Y(:,j))^10+100); % 可调
       % K2(i,j) = sqrt(norm(Y(:,i)-Y(:,j))^2+0.01); % 可调
       % K2(i,j) = 1- norm(Y(:,i)-Y(:,j))^2/(norm(Y(:,i)-Y(:,j))^2+0.1);
       % K2(i,j) = tanh(5*X(:,i)'*X(:,j)+1);
       K2(i,j) = (r+50*Y(:,i)'*Y(:,j))^d; % 可调
    end
end

for iter = 1 : 100
    H0 = H;
    H = C'*K;
    C0 = C;
    C1 = K*H';
    C2 = K*C*H*H';
    for i = 1 : n
        for j = 1 : k  
            C(i,j) = C(i,j)*C1(i,j)/C2(i,j);
        end
    end
    C = C./norm(C,'fro');
    norm(C-C0,'fro')/norm(C0,'fro');
    norm(H-H0,'fro')/norm(H0,'fro');
    iter;
    if max(norm(C-C0,'fro')/norm(C0,'fro'),norm(H-H0,'fro')/norm(H0,'fro')) < epslion
        iter
        break
    end
end

