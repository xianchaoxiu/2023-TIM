A = rand(2,2)

H = rand(2,2);
sigma=sqrt(1/2 * (H.^2*ones(2,2)))
abs(A(:))./sigma(:)