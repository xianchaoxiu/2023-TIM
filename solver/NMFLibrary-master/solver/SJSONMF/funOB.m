
%bsxfun+vecnorm
function B = funOB(B)
%global m k
m = 10;
k=5;
B = rand(m,k);
B = bsxfun(@rdivide,B,vecnorm(B));
end
