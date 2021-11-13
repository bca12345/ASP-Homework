%kappa = 0.667;
L = 1000;
R = 100;
x = zeros(100, 1000);
h = zeros(1, 1000);
f_1 = zeros(100, 1000);
S = load('ASP_HW2_Problem_4.mat');
P_f = zeros(1,10);
P_b = zeros(1,10);

% autocorrelation
r = zeros(1, 11);
for k = 0 : 10
    r(k+1) = -(128/105)*(1/4)^(abs(k)) + (64/21)*(1/2)^(abs(k));
end
[a, P, kappa] = ASP_Levison_Durbin(r);

for n = 0 : L-1
    h(n+1) = 2*(1/2)^n - (1/4)^n;
end
for k = 1 : R
    for i = 1 : L
        x(k,i) = fliplr(h(1:i))*S.V(k,1:i).'; % convolution
    end
end

%% forward prediction stage 1
for i = 1 : R
    f_1(i,:) = x(i,:) + kappa(2)*[x(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_1(n,:)).^2;
end
f1 = figure;
l = 1:1000;
plot(l,P_fn/R)
P_f(1) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 1
b_1 = zeros(100, 1000);

for i = 1 : R
    b_1(i,:) = [x(i,2:L) 0] + kappa(2)*x(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_1(n,:)).^2;
end
f2 = figure;
l = 1:1000;
plot(l,P_bn/R)
P_b(1) = sum(P_bn(250:750)/R) / 500;

%% forward prediction stage 2
f_2 = zeros(R,L);

for i = 1 : R
    f_2(i,:) = f_1(i,:) + kappa(3)*[b_1(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_2(n,:)).^2;
end
P_f(2) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 2
b_2 = zeros(R,L);

for i = 1 : R
    b_2(i,:) = [b_1(i,2:L) 0] + kappa(3)*f_1(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_2(n,:)).^2;
end
P_b(2) = sum(P_bn(250:750)/R) / 500;

%% forward prediction stage 3
f_3 = zeros(R,L);

for i = 1 : R
    f_3(i,:) = f_2(i,:) + kappa(4)*[b_2(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_3(n,:)).^2;
end
P_f(3) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 3
b_3 = zeros(R,L);

for i = 1 : R
    b_3(i,:) = [b_2(i,2:L) 0] + kappa(4)*f_2(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_2(n,:)).^2;
end
P_b(3) = sum(P_bn(250:750)/R) / 500;

%% forward prediction stage 4
f_4 = zeros(R,L);

for i = 1 : R
    f_4(i,:) = f_3(i,:) + kappa(5)*[b_3(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_4(n,:)).^2;
end
P_f(4) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 4
b_4 = zeros(R,L);

for i = 1 : R
    b_4(i,:) = [b_3(i,2:L) 0] + kappa(5)*f_3(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_4(n,:)).^2;
end
P_b(4) = sum(P_bn(250:750)/R) / 500;

%% forward prediction stage 5
f_5 = zeros(R,L);

for i = 1 : R
    f_5(i,:) = f_4(i,:) + kappa(6)*[b_4(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_5(n,:)).^2;
end
P_f(5) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 5
b_5 = zeros(R,L);

for i = 1 : R
    b_5(i,:) = [b_4(i,2:L) 0] + kappa(6)*f_4(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_5(n,:)).^2;
end
P_b(5) = sum(P_bn(250:750)/R) / 500;

%% forward prediction stage 6
f_6 = zeros(R,L);

for i = 1 : R
    f_6(i,:) = f_5(i,:) + kappa(7)*[b_5(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_6(n,:)).^2;
end
P_f(6) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 6
b_6 = zeros(R,L);

for i = 1 : R
    b_6(i,:) = [b_5(i,2:L) 0] + kappa(7)*f_5(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_6(n,:)).^2;
end
P_b(6) = sum(P_bn(250:750)/R) / 500;

%% forward prediction stage 7
f_7 = zeros(R,L);

for i = 1 : R
    f_7(i,:) = f_6(i,:) + kappa(8)*[b_6(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_7(n,:)).^2;
end
P_f(7) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 7
b_7 = zeros(R,L);

for i = 1 : R
    b_7(i,:) = [b_6(i,2:L) 0] + kappa(8)*f_6(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_7(n,:)).^2;
end
P_b(7) = sum(P_bn(250:750)/R) / 500;

%% forward prediction stage 8
f_8 = zeros(R,L);

for i = 1 : R
    f_8(i,:) = f_7(i,:) + kappa(9)*[b_7(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_8(n,:)).^2;
end
P_f(8) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 8
b_8 = zeros(R,L);

for i = 1 : R
    b_8(i,:) = [b_7(i,2:L) 0] + kappa(9)*f_7(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_8(n,:)).^2;
end
P_b(8) = sum(P_bn(250:750)/R) / 500;

%% forward prediction stage 9
f_9 = zeros(R,L);

for i = 1 : R
    f_9(i,:) = f_8(i,:) + kappa(10)*[b_8(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_9(n,:)).^2;
end
P_f(9) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 9
b_9 = zeros(R,L);

for i = 1 : R
    b_9(i,:) = [b_8(i,2:L) 0] + kappa(10)*f_8(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_9(n,:)).^2;
end
P_b(9) = sum(P_bn(250:750)/R) / 500;

%% forward prediction stage 10
f_10 = zeros(R,L);

for i = 1 : R
    f_10(i,:) = f_9(i,:) + kappa(11)*[b_9(i,2:L) 0]; % f(n) = x(n)+kappa*x(n-1)
end

P_fn = zeros(1,L);
for n = 1 : R
    P_fn = P_fn + abs(f_10(n,:)).^2;
end
P_f(10) = sum(P_fn(250:750)/R) / 500;

%% backward prediction stage 10
b_10 = zeros(R,L);

for i = 1 : R
    b_10(i,:) = [b_9(i,2:L) 0] + kappa(11)*f_9(i,:); % b(n) = x(n-1)+kappa*x(n)
end

P_bn = zeros(1,L);
for n = 1 : R
    P_bn = P_bn + abs(b_10(n,:)).^2;
end
P_b(10) = sum(P_bn(250:750)/R) / 500;

%% prediction error bound
S_x = @(f) log((-128/105)*( (1-(1/4)^2) ./ (1+(1/4)^2-2*(1/4)*cos(2*pi*f))) ...
    + (64/21)*( (1-(1/2)^2) ./ (1+(1/2)^2-2*(1/2)*cos(2*pi*f))));

q = integral(S_x,-(1/2),(1/2));
bound = exp(q);

axis = 1:10;
f3 = figure;
plot(axis,P_f,'-o',axis,P_b,'-x',axis,P(2:11),'-^',axis,bound*ones(1,length(axis)),'-*','LineWidth',2)
legend('Forward prediction error','Backward prediction error', ...
    'MMSE','Prediction error bound')
xlabel('Stage m')
ylabel('Prediction error')
title('ASP HW2 Problem 4f')