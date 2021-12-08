
% x is M-by-1 tap input vector
% d is desired response 
function [e,w] = ASP_LMS(x, d, mu)
M = 10;
w(M-1) = { zeros(M,1) };
e = zeros(1,500);

for n = M:length(x)
    x1 = x(n:-1:n-M+1);
    e(n) = d(n) - w{n-1}'*x1;
    w(n) = { w{n-1} + mu*conj(e(n))*x1 };
end