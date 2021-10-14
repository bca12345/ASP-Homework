%% the value of R, p, sd2, w1, w2, J_min
R = [2, 0.8, -0.4j; 0.8, 2, 0.8; 0.4j, 0.8, 2];
p = [1.6;-1.9;1.8];
sd2 = 12;
w1 = -0.5 + 1i; w2 = -1;
J_min = ASP_HW1_Wiener_MSE_5b(R, p, sd2);

%% find the relation between real(w0) and J(w)
w0_real = linspace(-4, 4, 201);
for j = 1 : 201
    w = [w0_real(j) + 1i; w1; w2];    
    J_w = ASP_Wiener_MSE(R, w, p, sd2);
    plot(w0_real(j), abs(J_w), 'b.-'); hold on
end
title('ASP_HW1_Wiener_MSE_5c');
xlabel('$Re\{w_0\}$','Interpreter','latex');
ylabel('MSE $|J|$', 'Interpreter','latex');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');


%% find extreme value
syms w0_real

w = [w0_real + 1i; w1; w2];
J_w0_real =  (w - R^(-1) * p)' * R * (w - R^(-1) * p) + J_min;
df = diff(J_w0_real, w0_real);
w0_min = solve(df==0, w0_real); % w0_min = 1
MSE = subs(J_w0_real, w0_real, w0_min);

text(w0_min, MSE, '\leftarrow extreme value')


