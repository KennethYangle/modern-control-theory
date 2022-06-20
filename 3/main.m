%% 鲁棒模型预测控制仿真
clear
close all
clc

%% 系统模型
A = [1.1, 2; 0, 0.95];
B = [0; 0.079];
C = eye(2);
D = [0.1; 0.5];

%% 控制器
% 控制时域
Nc = 10;

% 权重矩阵
Q = 1;
R = 0.1;
QQ = Q * eye(Nc * 2);
RR = R * eye(Nc * 1);

% 线性反馈
K = [2.5, 12.5];

% 求解P
Ak = A - B * K;
Qk = Q + K' * R * K;
P = dlyap(Ak', Qk);
QQ(Nc*2-1:Nc*2, Nc*2-1:Nc*2) = P;

% 计算状态序列矩阵
[F, Phi] = MPCGain(A, B, Nc);

% 约束
u_upper = 4;
u_lower = -4;

% 压缩约束
lb = (u_lower + 1.375) * ones(Nc, 1);
ub = (u_upper - 1.375) * ones(Nc, 1);

% 终端约束
u1_terminal_upper = 0.1;
u2_terminal_upper = 0.1;
u1_terminal_lower = -0.1;
u2_terminal_lower = -0.1;

% 仿真周期
total = 40;

% 优化变量
x = zeros(2, total+1);
u = zeros(1, total);
t = zeros(1, total);

% 标称系统变量
z = zeros(2, total+1);
v = zeros(1, total);

% 误差
e = zeros(2, total+1);

% 初始值
x(:, 1) = [1.1; -0.7];
z(:, 1) = x(:, 1);

%% 运行
for i = 1:total
    t(i) = i - 1;
    v0 = v(i);
    
    % 终端不等式约束
    bin1 =  [u1_terminal_upper; u2_terminal_upper] - A^Nc * z(:,1);
    bin2 = -[u1_terminal_lower; u2_terminal_lower] + A^Nc * z(:,1);
    bin = [bin1; bin2];
    Ain = [Phi(2*Nc-1:2*Nc, :); -Phi(2*Nc-1:2*Nc, :)];
    
    % 优化
    U = quadprog(Phi'*QQ*Phi+RR, z(:,i)'*F'*QQ*Phi, Ain, bin, [], [], lb, ub, v0);
    v(i) = U(1);
    e(:, i) = x(:, i) - z(:, i);
    u(i) = v(i) - K*e(:, i);
    
    % 模型迭代 
    w = (rand - 0.5) * 0.2;
    x(:, i+1) = A*x(:, i) + B*u(i) + D*w;
    
    % 标称模型
    z(:, i+1) = A*z(:, i) + B*v(i);
end

%% 画图
figure(1);
subplot(3, 1, 1);
plot(t, x(1, 1:total));
hold on;
grid on;

subplot(3, 1, 2);
plot(t, x(2, 1:total));
hold on;
grid on;

subplot(3, 1, 3);
plot(t, u);
hold off;