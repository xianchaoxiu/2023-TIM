function [ICn,Bn,d,ICc]=sort_IC(S,Q,B,cont)
W=B'*Q;
m=size(W,1);
for i=1:m
    w(i)=norm(W(i,:));
end;
[M,N]=sort(w,'descend');
for i=1:m
    ICa(i)=sum(M(1:i))/sum(M);
    if ICa(i)>=cont
        d=i;
        break;
    end;
end;
ICn=S(N(1:d),:);
Bn=B(:,N(1:d));
ICc=M/sum(M);