%% 鲁棒模型预测控制仿真
clear
close all
clc

%% 系统模型
A = [1, 1; 0, 1];
B = [0; 1];
C = eye(2);
D = [1; 1];

%% 控制器
% 控制时域
Nc = 20;

% 权重矩阵
Q = 1;
R = 0.1;
QQ = Q * eye(Nc * 2);
RR = R * eye(Nc * 1);

% 线性反馈
K = [0.4, 1.2];

% 求解P
Ak = A - B * K;
Qk = Q + K' * R * K;
P = dlyap(Ak', Qk);
QQ(Nc*2-1:Nc*2, Nc*2-1:Nc*2) = P;

% 计算状态序列矩阵
[F, Phi] = MPCGain(A, B, Nc);

% 约束
u_upper = 1;
u_lower = -1;

% 压缩约束
lb = (u_lower + 0.423) * ones(Nc, 1);
ub = (u_upper - 0.423) * ones(Nc, 1);

% 终端约束
u1_terminal_upper = 0.1;
u2_terminal_upper = 0.1;
u1_terminal_lower = -0.1;
u2_terminal_lower = -0.1;

% 仿真周期
total = 50;

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
X0 = [-1,1,1,0;-1,1,0,1];
x_stash = [];
u_stash = [];
for it = 1:4
    x(:, 1) = X0(:, it);
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
        w = (rand - 0.5) * 0.1;
        x(:, i+1) = A*x(:, i) + B*u(i) + D*w;

        % 标称模型
        z(:, i+1) = A*z(:, i) + B*v(i);
    end
    
    x_stash = [x_stash; x];
    u_stash = [u_stash; u];
end

%% 画图
figure(1);
subplot(3, 1, 1);
plot(t, x_stash(1, 1:total), t, x_stash(3, 1:total), t, x_stash(5, 1:total), t, x_stash(7, 1:total));
xlabel('迭代次数');
ylabel('x_1');
legend('[-1,-1]','[1,1]','[1,0]','[0,1]');
hold on;
grid on;

subplot(3, 1, 2);
plot(t, x(2, 1:total), t, x_stash(4, 1:total), t, x_stash(6, 1:total), t, x_stash(8, 1:total));
xlabel('迭代次数');
ylabel('x_2');
legend('[-1,-1]','[1,1]','[1,0]','[0,1]');
hold on;
grid on;

subplot(3, 1, 3);
plot(t, u_stash(1, 1:total), t, u_stash(2, 1:total), t, u_stash(3, 1:total), t, u_stash(4, 1:total));
xlabel('迭代次数');
ylabel('u');
legend('[-1,-1]','[1,1]','[1,0]','[0,1]');
hold off;
grid on;