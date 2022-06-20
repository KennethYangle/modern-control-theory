function [F, Phi] = MPCGain(A, B, N)
%MPC序列参数
[m1, n1] = size(A);
[n1, n_in] = size(B);

% 计算F
n = n1 + m1;
h(1:m1, :) = eye(m1);
F(1:m1, :) = eye(m1) * A;

for k = 2:N
    h(k*m1-1:k*m1, :) = h((k-1)*m1-1:(k-1)*m1, :) * A;
    F(k*m1-1:k*m1, :) = F((k-1)*m1-1:(k-1)*m1, :) * A;
end

% 计算Phi
v = h * B;
Phi = zeros(N*m1, N*n_in);
Phi(:, 1:n_in) = v;

for i = 2:N
    Phi(:, (i-1)*n_in+1:i*n_in) = [zeros((i-1)*m1, n_in); v(1:(N-i+1)*m1, 1:n_in)];
end

end

