function [I2,SPE]=variable_c(X,y,Devals,W,A)
n=size(y,2);
for i=1:n
    I2(i)=y(:,i)'*inv(Devals)*(y(:,i));%�������I2����ֵ
    SPE(i)=(X(:,i)-A*W*X(:,i))'*(X(:,i)-A*W*X(:,i));
end;