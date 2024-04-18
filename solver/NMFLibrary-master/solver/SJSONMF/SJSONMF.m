function [W,H] = SJSONMF(X)
global m n k s  
[m,n] = size(X);
s = 10;
W = rand(m,k);
Wf = S;
H = 10*rand(k, n);
tau1 = 10^-5;
tau2 = 10^-5;
V = ones(k,1)./sqrt(k);

%生成L
K = constructW0(X');
D = zeros(n);
for i = 1:n
    D(i,i) = sum(K(i,:));
end
L = D-K;

for iter = 1:500
    W0 = W;
    H0 = H;
    delta = 9.9*10^3; 
    tolf = 10^-4;
    epsg = 0.2;
    epsgm = 10^-4;
    gamma = 5;  
    eta = 0.8;  
    for witer = 1:500
        delta = 9.9*10^3;
        P = @(W) -1/delta*trace(X'*W*H)+norm(W*V,'fro')^2-1/delta*tau1*trace(W0'*W);
        sP = @(W) 2*W*V*V'-1/delta*X*H'-1/delta*tau1*W0;
        if P(W) > P(Wf)
            W = Wf;
        end
        Wbeta = 0.8;
        for sta = 1 : 5
            g1 = sP(W);
            d1 = -g1;
            PW = P(W);
            Walpha = 2;
            % Armijo line search 
            for v = 1:5
                intW = W+Walpha*d1;
                Wk = projOB(intW);
                tempW = norm(W - Wk,'fro')^2;
                Pk = P(Wk);
                Walpha = Walpha*Wbeta;   
                if Pk <= PW - Walpha/2*tempW 
                    break
                end
            end
            W = Wk;
            gP = sP(W)-W*diag(diag(W'*sP(W)));
            if norm((min(W,gP)),'fro') <= epsg && Pk <= PW %若Xt是否满足(4.1)和(4.2)则跳出循环 
                break
            end
        end 
        %终止条件
        norm(W*V,'fro')^2-1;
        if norm(W*V,'fro')^2-1 <= tolf
            break
        end
        delta = gamma*delta;
        epsg = max(eta*epsg,epsgm);
    end
    
    
    Hbeta = 0.5;
    Hepsilon = 10^-4;
    Hlambda = 10^-3;
    L2s = (1+tau2)/2+Hlambda*norm(L,'fro'); %需要计算
    fun = @(H) 0.5*norm(X-W*H,'fro')^2+Hlambda*trace(H*L*H')+0.5*tau2*norm(H-H0,'fro')^2;
    sfun = @(H) -W'*(X-W*H)+2*Hlambda*H*L+tau2*(H-H0);
    for hiter = 1 : 500
        f0 = fun(H);
        sf = sfun(H);
        Halpha = 1/L2s;  
        % Armijo line search 
        for q = 1 : 5
            intH = H - Halpha*sf;
            Hk = projH(intH);
            tempH = norm(H - Hk,'fro')^2;
            fk = fun(Hk);
            if fk <= f0 - Halpha/2*tempH
                break
            end
            Halpha = Halpha*Hbeta;  
        end
        H = Hk;
        d = find(H ~= 0); %给出非零位置
        [~,col] = find(H ~= 0); %给出非零位置行列，注意组合
        nabla_Gamma = sf(d); %给出全部非零元素
        if norm(nabla_Gamma,'fro') <= Hepsilon
            break; 
        end
    end
    werror = norm(W0-W,'fro')/norm(W0,'fro');
    herror = norm(H0-H,'fro')/norm(H0,'fro');
    if werror<=10^-5 && herror<=10^-5
       break
    end
end











