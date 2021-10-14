function [J_min] = ASP_HW1_Wiener_MSE_5b(R, p, sd2)
    w_opt = R^(-1) * p;
    J_min = sd2 - p' * w_opt;
end


