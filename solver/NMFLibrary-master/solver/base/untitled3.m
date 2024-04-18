global m n k 
X = rand(m,n);
[m,n]=size(X);
s = 11;
W0 = funOB;
Wf = S;
H = random('Normal', 1, 0.5, k, n);
delta = 10^-2; 
tolf = 10^-4;
epsg = 0.3;
epsgm = 0.2;
tmax = 30;
gamma2 = 5; 
eps = 0; 
gamma1 = 1; 
eta = 0.8;  
V = ones(k,1)./sqrt(k);
%P = @(W) 1/2*norm(X-W*H,'fro')^2+delta*(norm(W*V,'fro')^2-1);
% P = @(W) -1/delta*trace(X'*W*H)+norm(W*V,'fro')^2;
K = zeros(n,n);
for i = 1:n
    for j = 1:n
        K(i,j) = 1/(exp(norm(X(:,i)-X(:,j))^2)/0.3);
    end
end
D = zeros(n);
for i = 1:n
    D(i,i) = sum(sum(K(i,:)~=0));
end
L = D-K;
% E= -1/delta*trace(X'*W0*H)+norm(W0*V,'fro')^2+1/2*trace(H*L*H')
P = @(W) 0.5*norm(X-W*H,'fro')^2+delta*(norm(W*V,'fro')^2-1);
tic;
%外循环，若很多次数仍不满足 norm(Xt*V,'fro')^2-1 <= tolf，则跳出循环
for iter = 1:1000
    for witer = 1:tmax
        if P(W0) > P(Wf)
            W0 = Wf;
        else
            W0 = W0;
        end
        %调用梯度投影得到
        Agamma = 3;
        Abeta = 0.8;
        Asigma = 2;
        mmax = 8;
        %sP = @(W) -(X-W*H)*H'+2*delta*W*V*V';
        sP = @(W) 2*W*V*V'-1/delta*X*H';
        p = 30;
        PW0 = P(W0);
        for v = 1:mmax
            Wk = projOB(W0+Agamma*Abeta^v*d1);
            Pk = P(Wk);
            if Pk <= PW0 + Asigma*Agamma*Abeta^v*gd1
                break
            end
            if norm(Wk-W0,'fro') <= epsg && Pk <= PW0 %若Xt是否满足(4.1)和(4.2)则跳出循环            
                break
            end 
        end
        delta = gamma2*delta;
        if norm(Wk*V,'fro')^2-1 <= tolf
            break
        end
    end
    W0 = Wk;
    
    G = zeros(m,k);
    for i = 1:m
        w = max(W(i,:));
        j = find(W(i,:)==w);
        b = min(j);
        G(i,b) = 1;
    end
    WR = zeros(m,k);
    for j = 1:k
        WR(:,j) = (W(:,j).*G(:,j))/norm(W(:,j).*G(:,j));
    end
    WR(isnan(WR))=0;
    if vpa(WR'*WR) == eye(k,k)
        WR = WR;
    else
        WR = eye(m,k);
    end
    
    PCT = [];
    T = sign(WR);
    CT = X*H'.*T;
    for i = 1:k
        a = CT(:,i);
        e1 = zeros(m,1);
        e2 = zeros(m,1);
        b = find(a == max(a));%给出最大值点全部位置
        [z,c] = max(a);%给出第一个最大值点位置
        e1(sub2ind(size(e1), b)) = 1;
        e2(sub2ind(size(e2), c)) = 1;
        alpha1 = e1/norm(e1);
        alpha2 = e2/norm(e2);
        if z > 0
            a(a<0) = 0;
            Pai_a = a./norm(a);
        elseif z == 0
            Pai_a = alpha1;
        else
            Pai_a = alpha2;
        end
        PCT = [PCT,Pai_a];
    end  
    W0 = PCT;
    W1 = W0;

    M = 100;
    L2s = 1;%需要计算
    Hbeta = 0.8;
    Hepsilon = 10^-2;
    Hlamda = 2;
    fun = @(H) 0.5*norm(X-WR*H,'fro')^2;%+Hlamda*trace(H*L*H');
    sfun = @(H) -WR'*(X-WR*H);%+2*Hlamda*H*L;
    H00 = H;
    for hiter = 1 : M
        f0 = fun(H);
        sf = sfun(H);
        H0 = H;
        Halpha = 1/L2s;%给定几个值，看影响
        % Armijio line search 
        for q = 1 : 8
            intH = H0 - Halpha*sf;
            intH(intH<0) = 0;
            Z = intH;
            y = [];
            for i = 1:k
                yi = norm(Z(i,:));
                y = [y;yi];
            end
            [y, position] = maxk(y,s);
            H = zeros(size(Z));
            H(position,:) = Z(position,:);
            f = fun(H);
            temp = norm(H - H0,'fro')^2;
%             if f <= f0 + Hbeta/2*temp
%                 break; 
%             end
            Halpha = Halpha*Hbeta;  
        end
        f0 = f;
        d = find(H0 ~= 0 );%给出非零位置
        [row,col] = find(H0 ~= 0 );%给出非零位置行列，注意组合
        nabla_Gamma = sf(d);%给出全部非零元素
        if norm(nabla_Gamma) <= Hepsilon
            break; 
        end
    end

    eps = gamma1*eps;
    delta = gamma2*delta;
    epsg = max(eta*epsg,epsgm);
    time = toc;
end
W1 = W0;
H1 = H;



