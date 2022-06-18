function [u, p] = Controller(y, yd, dyd, ddyd, p)
%x state
%u control law

%% RBF
R = zeros(2, p.N);
for i = 1:p.N
    R(1,i) = 1/(sqrt(2*pi) * p.r) * exp( -(y-p.h(1,i))^2/(2*p.r^2) );
    R(2,i) = 1/(sqrt(2*pi) * p.r) * exp( -(y-p.h(2,i))^2/(2*p.r^2) );
end
R1 = R(1,:)';
R2 = R(2,:)';

H = [R1', zeros(1,p.N); zeros(1,p.N), R2'];
H1 = H(1,:);

dR = zeros(2, p.N);
for i = 1:p.N
    dR(1,i) = -(y-p.h(1,i))/p.r^2 * 1/(sqrt(2*pi) * p.r) * exp( -(y-p.h(1,i))^2/(2*p.r^2) );
    dR(2,i) = -(y-p.h(2,i))/p.r^2 * 1/(sqrt(2*pi) * p.r) * exp( -(y-p.h(2,i))^2/(2*p.r^2) );
end
dR1 = dR(1,:)';

%% 滤波器更新
K = [1; 1];
A0 = [-1,1; -1,0];

dxi = A0 * p.xi + K * y;
p.xi = p.xi + dxi * p.Ts;
xi2 = p.xi(2);
dXi = A0 * p.Xi + H;
p.Xi = p.Xi + dXi * p.Ts;
Xi2 = p.Xi(2,:);
dlambda = A0 * p.lambda + [0;1] * p.ModelParams.u;
p.lambda = p.lambda + dlambda * p.Ts;
lambda1 = p.lambda(1);
lambda2 = p.lambda(2);

%% alpha1及其各偏导数
% 第1步跟踪误差
z1 = y - yd;
% omega1
omega1 = [H1 + Xi2, 0]';
% alphabar
alpha1bar = -p.c1 * z1 - p.thetahat'*omega1 - xi2 - p.d1*z1 + dyd;
% p自适应律
dphat = p.gamma * z1 * alpha1bar - p.gamma * p.sigma2 * p.phat;
p.phat = p.phat + dphat * p.Ts;
% alpha1
alpha1 = p.phat * alpha1bar;
% 第2步跟踪误差
z2 = lambda2 - alpha1;
% alpha1的各偏导数
palpha1_pthetahat = -p.phat * omega1';
palpha1_py = -p.phat * ( p.c1 + p.d1 + p.thetahat' * [dR1; zeros(p.N+1,1) ] );
palpha1_pyd = p.phat * (p.c1 + p.d1);
palpha1_pdyd = p.phat;
%palpha1_pXi2 = -p.phat * p.thetahat' * [eye(2*p.N), zeros(2*p.N,1); zeros(1, 2*p.N+1)];
palpha1_pXi2 = -p.phat * p.thetahat(1:2*p.N,1);
palpha1_pphat = alpha1bar;

%% 调节函数和镇定函数
% 调节函数1
tau1 = p.Gamma * omega1 * z1 - p.sigma1 * p.Gamma * p.thetahat;
% omega2
omega2 = [palpha1_py * (Xi2 + H1), z1 - palpha1_py * lambda2]';
% 调节函数2
tau2 = tau1 + p.Gamma*omega2 * z2;
% beta2
beta2 = -lambda1 - palpha1_py * xi2 - (palpha1_pyd * dyd + palpha1_pdyd * ddyd + dXi(2,:) * palpha1_pXi2 + palpha1_pphat * dphat);

%% 自适应律和控制率
dthetahat = tau2;
p.thetahat = p.thetahat + dthetahat * p.Ts;
u = -p.c2*z2 - beta2 - p.thetahat'*omega2 - p.d2*z2*palpha1_py^2 + palpha1_pthetahat*tau2;
p.ModelParams.u = u;
end
