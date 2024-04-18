function SPE_control=control_limit(remain,c_aphla)
for i=1:3
    a(i)=sum(remain.^i);
end;
h=1-2*a(1)*a(3)/(3*a(2)^2);
SPE_control=a(1)*(c_aphla*h*sqrt(2*a(2))/a(1)+1+a(2)*h*(h-1)/a(1)^2)^(1/h);