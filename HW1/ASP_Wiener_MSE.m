function [J] = ASP_Wiener_MSE(R, w, p, sd2)
    %% check
    d = eig(R); %d(1) = lamda1, d(2) = lamda2, ...
    d_R = size(R);
    d_p = size(p);
    d_w = size(w);
    for i = 1 : 1 : d_R(1)
        if (d(i) < 0) || (d_R(1) ~= d_R(2))
            disp("error: Matrix is not positive semidefimite.")
        end
    end

    if (d_R(2) ~= d_w(1)) && (d_R(1) ~= d_p(1)) && (d_w(2) ~= d_p(2))
        disp("error: The dimensions of these input arguments are not suitable.")
    end

    if(imag(sd2) ~= 0)
        disp("error: The power is complex")
    end

    if(sd2 < 0)
        disp("error: The power is negative")
    end
    %% computation
    w_opt = R^(-1) * p;
    J_min = sd2 - p' * w_opt;
    J = (w - R^(-1) * p)' * R * (w - R^(-1) * p) + J_min;
end
