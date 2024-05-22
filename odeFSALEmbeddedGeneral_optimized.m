function [t, xsol, ELE] = odeFSALEmbeddedGeneral_optimized(c_vector, A_matrix, b_vector, b_hat_vector, odefun, tspan, tau, incond)

s_stages = length(c_vector);
m = length(incond);

c_vector = reshape(c_vector, [s_stages 1]);
b_vector = reshape(b_vector, [s_stages 1]);
b_hat_vector = reshape(b_hat_vector, [s_stages 1]);
d_vector = b_hat_vector - b_vector;
incond = reshape(incond, [m 1]);

t = (tspan(1) : tau : tspan(2))';
xsol = zeros(length(incond), length(t));
xsol(:, 1) = incond(:);
K_matrix = zeros(m, s_stages);
ELE = zeros(length(t), 1);
K_matrix(:, s_stages) = odefun(t(1), xsol(:, 1));

for n = 1:length(t)-1
    K_matrix(:, 1) = K_matrix(:, s_stages);
        for i = 2:s_stages
            K_matrix(:, i) = odefun(t(n) + tau * c_vector(i), xsol(:, n) + tau * K_matrix(:, 1:i-1) * A_matrix(i, 1:i-1)');
        end
    xsol(:, n+1) = xsol(:, n) + tau * K_matrix * b_vector;
    ELE(n+1) = norm(tau * K_matrix * d_vector, "inf");
end
xsol = xsol';
end