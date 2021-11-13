function[a, P, kappa] = ASP_Levison_Durbin(r)
   
    P(1) = r(1);
    c(1) = conj(r(2));
    a = {1}; % a_0 = 1
    R = {r(2)};
    for i = 2 : length(r)
        kappa(i) = -(c(i - 1) / P(i - 1));
        a(i) = {[a{i-1} 0] + [0 fliplr(kappa(i)*a{i-1})]};
        P(i) = P(i - 1) * (1 - abs(kappa(i))^2);
        if i+1 < length(r) + 1
            R(i) = {[r(i+1) R{i-1}]};
            c(i) = R{i} * a{i}';
        end
    end
end
