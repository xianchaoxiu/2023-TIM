function [yn,Bn,d,Devals,ICc]=sort_IC_cum(y,B,evals,cont)
Evals=diag(evals);
m=length(Evals);
[M,N]=sort(Evals,'descend');
for i=1:m
    ICa(i)=sum(M(1:i))/sum(M);
    if ICa(i)>=cont
        d=i;
        break;
    end;
end;
Devals=diag(Evals(N(1:d)));
yn=y(N(1:d),:);
Bn=B(:,N(1:d));
ICc=M/sum(M);