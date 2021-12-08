%% problem 4(c)
M = 10;
sd2 = 2;
mu = 0.02;
% impulse response
% h = zeros(1, 500);
% for n = 0 : L-1
%     h(n+1) = (5/6)*(1/2)^n + (1/6)*(-1/4)^n;
% end
% x = zeros(R,L);
% for k = 1 : R
%     for i = 1 : L
%         x(k,i) = fliplr(h(1:i))*S.matV(k,1:i).'; % convolution
%     end
% end

% compute optimal weight vector in Wiener filters
% for i = 0 : 9
%     p(i+1) = (5/6)*(1/2)^i + (1/6)*(-1/4)^i;
% end
% p = [h(0);h(-1);...;h(-M+1)] and h(n) = 0 for n < 0
p = [1; zeros(M-1,1)];

% autocorrealation matrix
fun = @(f,n)( (1 + 1/64 + (1/8) * (exp(1j*2*pi*f)+exp(-1j*2*pi*f))) / ...
    (1 + 1/16 + 1/64 + (-7/32) * (exp(1j*2*pi*f)+exp(-1j*2*pi*f)) + ...
    (1/8) * (exp(1j*4*pi*f) + exp(1j*4*pi*f)) ) ) * exp(1j*2*pi*f*n);
F =@(n)integral(@(f)fun(f,n),0,1);

next = 0;
r = zeros(M,M);
for i = 0 : M-1
    for j = 0 : M-1
        r(i+1,j+1) = F(j - next);
    end
    next = next + 1;
end
[U,V] = eig(r); % V: eigenvalues, U: eigenvectors
eig_value = diag(V);

% optimal weight vector
w_opt = r^(-1) * p;

% steepest descent
% nu = cell(500);
w = cell(500);
w(1) = { zeros(M,1) };
% c = cell(500);
nu_0 = U'*w_opt;

for n = 1 : 500
%     c(n) = { w_opt - w{1} }; % weight error vector
%     nu(n) = { U'*c{1} };     % nu(n) = [nu_1(n); ... ;nu_M(n)]

    % MSE
    J_min = sd2 - p' * w_opt;
    s = 0;
    for k = 1:M
        s = s + eig_value(k)*(1-mu*eig_value(k))^(2*n)*abs(nu_0(k))^2;
    end
    J(n) = J_min + s;

    % MSD    
    if n < 500
    e_msd = w_opt - w{n};
    w(n+1) = { (eye(M) - mu*r)*w{n} + mu*p };
    D(n) = norm(e_msd)^2;
    end
end
plot([1:500],abs(J)); hold on
plot([1:499],abs(D));
title('ASP HW3 Problem 4c')
xlabel('Number of adaption cycle n');
ylabel('error');
legend('MSE learning curve', 'MSD learning curve');
hold off;
% function [J_min] = Jmin(R, p, sd2)
%     w_opt = R^(-1) * p;
%     J_min = sd2 - p' * w_opt;
% end
