function [e,w] =  ASP_RLS(x, d, lambda, delta)
M = 10;
w(M-1) = { zeros(M,1) }; % use cell array
e = zeros(1,500);
P = (1/delta) * eye(M); % P matrix

for n = M:length(x)
    x1 = x(n:-1:n-M+1); % note x1 is 10x1
    k = (lambda^(-1)*P*x1) / (1 + lambda^(-1)*x1'*P*x1);
    e(n) = d(n) - w{n-1}'*x1;
    w(n) = { w{n-1} + k * conj(e(n)) };
    P = (P - k*x1'* P) / lambda;
end