%% 仿真设置
p.simN = 5000;            % 迭代次数
p.Ts = 0.01;              % 迭代步长
p.plotOn = true;          % 画图开关

%% 系统参数
p.ModelParams.b = -1;     % 系统参数
p.ModelParams.x0 = [1;0];
p.ModelParams.u = 0;

%% RBF参数
p.r = 2;
p.h = [-4,-3,-2,-1,0,1,2,3,4; -4,-3,-2,-1,0,1,2,3,4];
[~, p.N] = size(p.h);
p.thetahat = zeros(2*p.N+1, 1);

%% 滤波器初值
p.xi = zeros(2,1);
p.Xi = zeros(2,2*p.N);
p.lambda = zeros(2,1);
p.phat = -1;

%% 设计参数
p.c1 = 2;
p.d1 = 2;
p.c2 = 1;
p.Gamma = 3*eye(2*p.N+1);
p.gamma = 3;
p.sigma1 = 1;
p.sigma2 = 1;