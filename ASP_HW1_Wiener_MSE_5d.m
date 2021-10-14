%% the value of R, p, sd2, w1, w2, J_min
R = [2, 0.8, -0.4j; 0.8, 2, 0.8; 0.4j, 0.8, 2];
p = [1.6;-1.9;1.8];
sd2 = 12;
w1 = -2.5 + 0.1i; w2 = 2 - 0.4i;
J_min = ASP_HW1_Wiener_MSE_5b(R, p, sd2);
J_e_w = zeros(201, 201);
J_w = zeros(201, 201);

%% quadratic form
w0_real = linspace(-4, 4, 201);
w0_imag = linspace(-4, 4, 201);

for m = 1 : 201
    for n = 1 : 201 
    w = [w0_real(m) + w0_imag(n)*1i; w1; w2];
    J_w(m, n) = ASP_Wiener_MSE(R, w, p, sd2);
    end
end

%% figure

title('ASP HW1 Problem 5c');
m = 1:201; n = 1:201;
[M, N] = meshgrid(m, n);
h = surf(M, N, abs(J_w));
xlabel('$Re\{w_0\} = -4 + 0.04X$','Interpreter','latex');
ylabel('$Im\{w_0\} = -4 + 0.04Y$', 'Interpreter','latex');
zlabel('MSE $|J|$', 'Interpreter','latex');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
set(h, 'LineStyle', 'none')
colorbar

%% find extreme value

syms w0_real real
syms w0_imag 

w = [w0_real + w0_imag*1i ; w1 ; w2];
J_w0 =  ASP_Wiener_MSE(R, w, p, sd2);
df1 = diff(J_w0, w0_real);
df2 = diff(J_w0, w0_imag);
df = [0.5*df1 0.5*df2];
s = solve(df, [w0_real w0_imag]);
MSE = subs(J_w0, [w0_real w0_imag], [s.w0_real s.w0_imag]);

text((s.w0_imag + 4)/0.04, (s.w0_real + 4)/0.04, MSE, '\leftarrow extreme value')

