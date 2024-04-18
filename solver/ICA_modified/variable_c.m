function [I2,SPE]=variable_c(X,y,Devals,W,A)
n=size(y,2);
for i=1:n
    I2(i)=y(:,i)'*inv(Devals)*(y(:,i));%π ’œ’Ô∂œI2π±œ◊÷µ
    SPE(i)=(X(:,i)-A*W*X(:,i))'*(X(:,i)-A*W*X(:,i));
end;