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
end

%% 画图
if p.plotOn == true
    t = p.Ts * (1: p.simN);
    figure(1);
    plot(t, x_stash(1,:)/1.5, t, x_stash(2,:)/1.5);
    legend('x_1','x_2');
    figure(2);
    plot(t, x_stash(1,:)/1.5, t, yd_stash(:));
    legend('y','y_d');
end
