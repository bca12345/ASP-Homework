S = load("ASP_Final_Data.mat");

A = S.matX;
[r, c] = size(A);
theta = -90:1:90;
theta_ss = zeros(1, length(A));
theta_ii = zeros(1, length(A));
s_t_hat = zeros(1, length(A)); % source_lms
source_mvdr = zeros(1, length(A));
d = zeros(1, length(A));

w_mvdr = zeros(r, c);
w_lms = zeros(r, c);
s = 0;

% estimate theta
for i = 20 : 20 : length(A)
    [phi1, phi1_inv] = correlation_matrix_RLS(A(:, 1+s:i), 0.9);
    P = mvdr_spectrum(phi1, phi1_inv);
    [v, l] = findpeaks(10*log10(abs(P)), 'MinPeakHeight', abs(mean(10*log10(abs(P))))*0.95 );
    d(i) = length(v);
    s = i;
    if length(v) == 2
        theta_ss(i) = max(l) - 91;
        theta_ii(i) = min(l) - 91;
    elseif length(v) == 1    %only one peak
        if l(1) > 90 && l(1) < 100
            theta_ss(i) = l - 91;
            theta_ii(i) = -l + 91;
        else
            theta_ii(i) = l - 91;
            theta_ss(i) = -l + 91;
        end
    end
end
theta_s_hat = replace_zero(theta_ss);
theta_i_hat = replace_zero(theta_ii);
% plot(theta, 10*log10(abs(P)))

% w_mvdr and w_lms
s = 0;
pos = 0 : 9;
for i = 20 : length(A)
    [phi1, phi1_inv] = correlation_matrix_RLS(A(:, 1+s:i), 0.9);
    
    a = exp(1j*pi*sind(theta_s_hat(i).*pos(:)));
    w_mvdr(:, i) = phi1_inv*a/(a'*phi1_inv*a);
    w_lms(:, i) = lms_weight(A(:, 1+s:i), 1 * theta_s_hat(i), 1 * theta_i_hat(i));
    s = s + 1;
    s_t_hat(i) = w_lms(:, i)' * A(:, i);
    source_mvdr(i) = w_mvdr(:, i)' * A(:, i);
end

%w_uni
w_uni = zeros(r, c);
source_uni = zeros(1, c);
pos = 0 : 9;
for i = 1 : length(A)
    w_uni(:, i) = (1/r) * exp(1j*pi*sind(theta_s_hat(i).*pos(:)));
    source_uni(i) = w_uni(:, i)' * A(:, i);
end

% w_lcmv
w_lcmv = zeros(10, length(A));
source_lcmv = zeros(1, length(A));
pos = 0:9;
s = 0;
for i = 20 : length(A)
    [phi1, phi1_inv] = correlation_matrix_RLS(A(:, 1+s:i), 0.9);
    sv(:, 1) = exp(1j*pi*sind(theta_s_hat(i).*pos(:)));
    sv(:, 2) = exp(1j*pi*sind(theta_i_hat(i).*pos(:)));
    resp = [1; 0];
    w_lcmv(:, i) = phi1_inv * sv * (sv'*phi1_inv*sv + 1e-4*eye(2))^(-1) * resp;
    s = s + 1;
    source_lcmv(i) = w_lcmv(:, i)' * A(:, i);
end

figure
plot(abs(source_mvdr))
title('MVDR')
xlabel('time index')
ylabel('$|\hat{s}(t)|$','interpreter','latex')

figure
plot(abs(s_t_hat))
title('LMS with linear constraint')
xlabel('time index')
ylabel('$|\hat{s}(t)|$','interpreter','latex')

figure
plot(theta_s_hat)
title('the source DOA')
xlabel('time index')
ylabel('$\hat{\theta}_s$','interpreter','latex')

figure
plot(theta_i_hat)
title('the interference DOA')
xlabel('time index')
ylabel('$\hat{\theta}_i$','interpreter','latex')

save('theta_s_hat.mat', "theta_s_hat")
save('theta_i_hat.mat', "theta_i_hat")
save('s_t_hat.mat', "s_t_hat")

% function
function [phi, phi_inv] = correlation_matrix_RLS(input_matrix, lambda)
    [r, c] = size(input_matrix);
    phi = zeros(r, r);
    phi_inv = (1/0.01) * eye(r);
    for m = 1 : length(input_matrix)
        x = input_matrix(:, m);
        k = (lambda^(-1)*phi_inv*x) / (1 + lambda^(-1)*x'*phi_inv*x);
        phi_inv = (phi_inv - k*x'* phi_inv) / lambda;
        phi = lambda * phi + (x * x') + 1e-4 * eye(r);
    end
end

function [P] = mvdr_spectrum(input_matrix, inverse_matrix)
    theta = -90:1:90;
    N = length(input_matrix);
    P = zeros(1, length(theta));
    a = zeros(N, 1);
    position = 0:9;
    for i = 1 : length(theta)
        for n = 1 : N
            a(n, 1) = exp(1j*pi*sind(theta(i))*position(n));
        end
%         w_mvdr = inverse_matrix*a/(a'*inverse_matrix*a);
%         ww(i, :) = input_matrix^(-1)*a/(a'*input_matrix^(-1)*a);
        P(i) = 1/(a'*inverse_matrix*a);
    end
end

function [theta_hat] = replace_zero(input_theta)
    for i = 2 : length(input_theta)
        if input_theta(i) == 0
            input_theta(i) = input_theta(i - 1);
        end
    end
    theta_hat = input_theta;
end

function [w_lms] = lms_weight(input_matrix, theta_s, theta_i)
%   A = S.matX;
%     if theta_s == theta_i
%         theta_i = theta_i - 4;
%     end
    [r, c] = size(input_matrix);
    position = 0:9;
    C(:, 1) = exp(1j*pi*sind(theta_s).*position(:)); % constraint matrix
    C(:, 2) = exp(1j*pi*sind(theta_i).*position(:));
    C_a = null(C'); % orthogonal complement    
    g = [1; 1e-4]; % gain vector
    w_q = C * (C' * C)^(-1) * g;
    M = 10;
    w_a = zeros(r - length(C(1, :)), 1);
    mu = 0.001;
    w_lms = zeros(r, 1);
    for n = M:length(input_matrix)
        x = input_matrix(:, n:-1:n-M+1);
        w_a = w_a + mu * C_a' * (x * x') * (w_q - C_a * w_a);
    %     P(n) = w_mvdr'*input_matrix*w_mvdr;
        w_lms = w_q - C_a * w_a;
    end
end
