% Copyright (C) 2018, ETH Zurich, D-ITET, Kenneth Kuchera, Alexander Liniger
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xp=SimTimeStep(x,u,p)
%x state
%u input
%Ts sampling time
x0=x;
[~, inivt]=ode45(@(t, x)Sys(t, x, u, p.ModelParams), [0, p.Ts], x0);
xp=inivt(end,:);
end

function f1y=f1(y)
f1y = y^2;
end

function f2y=f2(y)
f2y = 2*y;
end

% 系统方程
function xdot=Sys(t, x, u, ModelParams)
b = ModelParams.b;
y = x(1);

xdot = [x(2) + f1(y);
        b*u + f2(y)];
end
