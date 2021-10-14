%% the value of R, p, sd2, w1, w2, J_min
R = [2, 0.8, -0.4j; 0.8, 2, 0.8; 0.4j, 0.8, 2];
p = [1.6;-1.9;1.8];
sd2 = 12;
w2 = 2 + -0.37*1i;
J_min = ASP_HW1_Wiener_MSE_5b(R, p, sd2);
J_e_w = zeros(201, 201);
J_w = zeros(201, 201);
[U, D] = eig(R);

%% compute quadratic form
w0_r = linspace(-3, 3, 201);
w1_r = linspace(-3, 3, 201);

for m = 1 : 201
    for n = 1 : 201
        w = [w0_r(m) + 0.4*1i ; w1_r(n) - 0.00125*1i ; w2];
        J_w(m, n) = ASP_Wiener_MSE(R, w, p, sd2);
    end
end

%% figure
[M, N] = meshgrid(1:201, 1:201);
lev = [0.8, 0.9, 1, 2, 3, 4, 5];
contour(M, N, abs(J_w), lev, 'ShowText','on')
xlabel('$Re\{w_0\} = -3 + 0.03X$','Interpreter','latex');
ylabel('$Re\{w_1\} = -3 + 0.03Y$', 'Interpreter','latex');

%% find extreme value
syms w0_real real
syms w1_real real

w = [w0_real + 0.4*1i ; w1_real - 0.00125*1i ; w2];
J_w0_w1 =  ASP_Wiener_MSE(R, w, p, sd2);

df1 = diff(J_w0_w1, w0_real) ;
df2 = diff(J_w0_w1, w1_real);
df = [0.5*df1 0.5*df2];
s = solve(df, [w0_real w1_real]);

text( (s.w1_real + 3)/0.03, (s.w0_real + 3)/0.03, '\leftarrow extreme value')

