%% RBF鲁棒自适应控制仿真
clear
close all
clc

%% 载入参数
LoadParams;

%% 仿真
x = p.ModelParams.x0;
x_stash = zeros(2, p.simN);
yd_stash = zeros(1, p.simN);
u_stash = zeros(1, p.simN);
for i = 1: p.simN
    % 控制器
    y = x(1);
    [yd, dyd, ddyd] = Reference(i * p.Ts);
    [u, p] = Controller(y, yd, dyd, ddyd, p);
    % 系统运行
    x = SimTimeStep(x, u, p)';
    % 记录
    x_stash(:, i) = x;
    yd_stash(i) = yd;
    u_stash(i) = u;
end

%% 画图
if p.plotOn == true
    t = p.Ts * (1: p.simN);
    figure(1);
    subplot(3, 1, 1);
    plot(t, x_stash(1,:)/1.5, t, yd_stash(:));
    legend('y','y_d');
    xlabel('t / s');
    ylabel('y');
    grid on;
    
    subplot(3, 1, 2);
    plot(t, x_stash(1,:)/1.5, t, x_stash(2,:)/1.5);
    xlabel('t / s');
    ylabel('x');
    legend('x_1','x_2');
    grid on;
    
    subplot(3, 1, 3);
    plot(t, u_stash);
    xlabel('t / s');
    ylabel('u');
    legend('u');
    grid on;
end
