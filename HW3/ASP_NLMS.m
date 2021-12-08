function [e,w] = ASP_NLMS(x, d, mu)
M = 10;
w(M-1) = { zeros(M,1) };
e = zeros(1,500);

for n = M:length(x)
    x1 = x(n:-1:n-M+1);
    e(n) = d(n) - w{n-1}'*x1;
    nn = norm(x1)^(-2);
    w(n) = { w{n-1} + mu*conj(e(n))*x1*nn };
end