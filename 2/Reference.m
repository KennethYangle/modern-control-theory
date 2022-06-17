function [yd, dyd, ddyd]=Reference(t)
% 参考信号
yd = sin(0.5*t);
dyd = 0.5*cos(0.5*t);
ddyd = -0.25*sin(0.5*t);
%yd = 0.5*t*3;
%dyd = 1.5*t^2;
%ddyd = 3*t;
end