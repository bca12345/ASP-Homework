%% problem 4(b)
L = 500;
R = 1000;
S = load('ASP_HW3_Problem_4.mat');
M = 10;
% impulse response
h = zeros(1, 500);
for n = 0 : L-1
    h(n+1) = (5/6)*(1/2)^n + (1/6)*(-1/4)^n;
end
x = zeros(R,L);
for k = 1 : R
    for i = 1 : L
        x(k,i) = fliplr(h(1:i))*S.matV(k,1:i).'; % convolution
    end
end

% compute optimal weight vector in Wiener filters
% p = [h(0);h(-1);...;h(-M+1)] and h(n) = 0 for n < 0
p = [1; zeros(M-1,1)];

% fun_h = @(f,n)( (1 + (1/8)*exp(-1j*2*pi*f))/(1 - (1/4) * exp(-1j*2*pi*f) - (1/8) * exp(-1j*4*pi*f)) )...
%     * exp(1j*2*pi*f*n);
% H =@(n)integral(@(f)fun(f,n),-1/2,1/2);
% 
% p = zeros(M,1);
% for i = 0 : M-1
%     p(i+1) = conj(H(i));
% end

% autocorrealation matrix
fun = @(f,n)( (1 + 1/64 + (1/8) * (exp(1j*2*pi*f)+exp(-1j*2*pi*f))) / ...
    (1 + 1/16 + 1/64 + (-7/32) * (exp(1j*2*pi*f)+exp(-1j*2*pi*f)) + ...
    (1/8) * (exp(1j*4*pi*f) + exp(1j*4*pi*f)) ) ) * exp(1j*2*pi*f*n);
F =@(n)integral(@(f)fun(f,n),-1/2,1/2);

next = 0;
r = zeros(M,M);
for i = 0 : M-1
    for j = 0 : M-1
        r(i+1,j+1) = F(j - next);
    end
    next = next + 1;
end

% optimal weight vector
w_opt = r^(-1) * p;

% LMS with mu = 0.02
d = S.matV(1:R,:);
wn = cell(R,L);
e = zeros(1,L);
for j = 1 : R
    mu = 0.02;
    [e(j,:), wn(j,:)] = ASP_LMS(x(j,:).',d(j,:),mu);
end

% MSD
e_msd = cell(R,L);
for i = 1 : R
    for j = M : L
        e_msd(i,j) = { w_opt - wn{i,j} };
    end
end

for i = 1 : L
    s = 0;
    for j = 1 : R
        s = s + norm(e_msd{j,i})^2;
    end
    D1(i) = (1/R) * s;
end

% MSE
for i = 1 : L
    J1(i) = (1/R) * sum(abs(e(:,i)).^2);
end

% LMS with mu = 0.1
for j = 1 : R
    mu = 0.1;
    [e(j,:), wn(j,:)] = ASP_LMS(x(j,:).',d(j,:),mu);
end

% MSD
e_msd = cell(R,L);
for i = 1 : R
    for j = M : L
        e_msd(i,j) = { w_opt - wn{i,j} };
    end
end

for i = 1 : L
    s = 0;
    for j = 1 : R
        s = s + norm(e_msd{j,i})^2;
    end
    D2(i) = (1/R) * s;
end

% MSE
for i = 1 : L
    J2(i) = (1/R) * sum(abs(e(:,i)).^2);
end

% NLMS with mu = 0.3
for j = 1 : R
    mu = 0.3;
    [e(j,:), wn(j,:)] = ASP_NLMS(x(j,:).',d(j,:),mu);
end

% MSD
e_msd = cell(R,L);
for i = 1 : R
    for j = M : L
        e_msd(i,j) = { w_opt - wn{i,j} };
    end
end

for i = 1 : L
    s = 0;
    for j = 1 : R
        s = s + norm(e_msd{j,i})^2;
    end
    D3(i) = (1/R) * s;
end

% MSE
for i = 1 : L
    J3(i) = (1/R) * sum(abs(e(:,i)).^2);
end

% NLMS with mu = 0.9
for j = 1 : R
    mu = 0.9;
    [e(j,:), wn(j,:)] = ASP_NLMS(x(j,:).',d(j,:),mu);
end

% MSD
e_msd = cell(R,L);
for i = 1 : R
    for j = M : L
        e_msd(i,j) = { w_opt - wn{i,j} };
    end
end

for i = 1 : L
    s = 0;
    for j = 1 : R
        s = s + norm(e_msd{j,i})^2;
    end
    D4(i) = (1/R) * s;
end

% MSE
for i = 1 : L
    J4(i) = (1/R) * sum(abs(e(:,i)).^2);
end

% RLS with lambda = 0.8, delta = 0.01
lambda = 0.8; delta = 0.01;
for j = 1 : R
    [e(j,:), wn(j,:)] = ASP_RLS(x(j,:).',d(j,:),lambda,delta);
end

% MSD
e_msd = cell(R,L);
for i = 1 : R
    for j = M : L
        e_msd(i,j) = { w_opt - wn{i,j} };
    end
end

for i = 1 : L
    s = 0;
    for j = 1 : R
        s = s + norm(e_msd{j,i})^2;
    end
    D5(i) = (1/R) * s;
end

% MSE
for i = 1 : L
    J5(i) = (1/R) * sum(abs(e(:,i)).^2);
end

% RLS with lambda = 0.8, delta = 0.1
lambda = 0.8; delta = 0.1;
%wn = cell(R,L);
for j = 1 : R
    [e(j,:), wn(j,:)] = ASP_RLS(x(j,:).',d(j,:),lambda,delta);
end

% MSD
e_msd = cell(R,L);
for i = 1 : R
    for j = M : L
        e_msd(i,j) = { w_opt - wn{i,j} };
    end
end

for i = 1 : L
    s = 0;
    for j = 1 : R
        s = s + norm(e_msd{j,i})^2;
    end
    D6(i) = (1/R) * s;
end

% MSE
for i = 1 : L
    J6(i) = (1/R) * sum(abs(e(:,i)).^2);
end


% RLS with lambda = 0.99, delta = 0.01
lambda = 0.99; delta = 0.01;
for j = 1 : R
    [e(j,:), wn(j,:)] = ASP_RLS(x(j,:).',d(j,:),lambda,delta);
end

% MSD
e_msd = cell(R,L);
for i = 1 : R
    for j = M : L
        e_msd(i,j) = { w_opt - wn{i,j} };
    end
end

for i = 1 : L
    s = 0;
    for j = 1 : R
        s = s + norm(e_msd{j,i})^2;
    end
    D7(i) = (1/R) * s;
end

% MSE
for i = 1 : L
    J7(i) = (1/R) * sum(abs(e(:,i)).^2);
end

semilogy(J1,'b'); hold on
semilogy(J2,'k');
semilogy(J3,'r');
semilogy(J4,'g');
semilogy(J5,'.-');
semilogy(J6,'--');
semilogy(J7,':');
legend('LMS with \mu = 0.02', ...
    'LMS with \mu = 0.1', ...
    'NLMS with \mu = 0.3', ...
    'NLMS with \mu = 0.9', ...
    'RLS with \lambda = 0.8, \delta = 0.01', ...
    'RLS with \lambda = 0.8, \delta = 0.1', ...
    'RLS with \lambda = 0.99, \delta = 0.01')

title('MSE learning curves with R=1000')
xlabel('Number of adaption cycles n')
ylabel('Mean-square error')
hold off;

figure
semilogy(D1,'b'); hold on
semilogy(D2,'k');
semilogy(D3,'r');
semilogy(D4,'g');
semilogy(D5,'.-');
semilogy(D6,'--');
semilogy(D7,':');
legend('LMS with \mu = 0.02', ...
    'LMS with \mu = 0.1', ...
    'NLMS with \mu = 0.3', ...
    'NLMS with \mu = 0.9', ...
    'RLS with \lambda = 0.8, \delta = 0.01', ...
    'RLS with \lambda = 0.8, \delta = 0.1', ...
    'RLS with \lambda = 0.99, \delta = 0.01')

title('MSD learning curves with R=1000')
xlabel('Number of adaption cycles n')
ylabel('Mean-square deviation')
hold off;
