clc;clear all;close all;

%常量定义
m = 1.6726e-27;
q = 1.602e-19;
eV = 1.602e-19;
velocity = sqrt((2 * 10000 * eV) / m);
v0 = [velocity * cos(pi / 4); velocity * cos(pi / 3); velocity * sin(pi / 3)];
r0 = [5 * (sqrt(2) / 2) * 6.371e6; 5 * (sqrt(2) / 2) * 6.371e6; 0];
Vacuum_permeability = 4 * pi * 1e-7;
dipole_moment = 7.94e22;
B = get_magneticfield(Vacuum_permeability , dipole_moment , r0);


%步长周期
T = 6e4;
dt = 5e-4;
steps = T / dt;

%数据记录
trace = zeros(3, steps+1);
trace(:,1) = r0;

%RK4
r = r0;
v = v0;
for step = 1:steps
    %K1
    [a1, dr1] = lorentz_force(r, v, q, m, B);
    
    % ---- k2 计算 ----
    v_temp = v + 0.5*dt*a1;
    r_temp = r + 0.5*dt*dr1;
    [a2, dr2] = lorentz_force(r_temp, v_temp, q, m, B);
    
    % ---- k3 计算 ----
    v_temp = v + 0.5*dt*a2;
    r_temp = r + 0.5*dt*dr2;
    [a3, dr3] = lorentz_force(r_temp, v_temp, q, m, B);
    
    % ---- k4 计算 ----
    v_temp = v + dt*a3;
    r_temp = r + dt*dr3;
    [a4, dr4] = lorentz_force(r_temp, v_temp, q, m, B);
    
    % ---- 更新 ----
    v = v + dt/6*(a1 + 2*a2 + 2*a3 + a4);
    r = r + dt/6*(dr1 + 2*dr2 + 2*dr3 + dr4);
    
    trace(:,step+1) = r;
end

function [accel, velocity] = lorentz_force(r, v, q, m, B)
    force = q * cross(v, B);
    accel = force/m;    % acc
    velocity = v;       % vel
end

function B = get_magneticfield(Vacuum_permeability , dipole_moment , r)
    B(1,:) = ((3 * Vacuum_permeability * dipole_moment) / (4 * pi * norm(r) ^5 )) * r(1) * r(3);
    B(2,:) = ((3 * Vacuum_permeability * dipole_moment) / (4 * pi * norm(r) ^5 )) * r(2) * r(3);
    B(3,:) = ((3 * Vacuum_permeability * dipole_moment) / (4 * pi * norm(r) ^5 )) * (3 * r(3) ^2 - norm(r) ^2);
end
%绘图
x = trace(1,:);
y = trace(2,:);
z = trace(3,:);

% 投影到 XY 平面
figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1);  
plot(x, y, "r-", "LineWidth", 1);
xlabel('X');
ylabel('Y');
title('XY Projection');
axis equal;
grid on;

% 投影到 XZ 平面
subplot(1,3,2);
plot(x, z, "g-", "LineWidth", 1);
xlabel('X');
ylabel('Z');
title('XZ Projection');
grid on;

% 投影到 YZ 平面
subplot(1,3,3);
plot(y, z, "b-", "LineWidth", 1);
xlabel('Y');
ylabel('Z');
title('YZ Projection');
grid on;

%粒子全局磁场
figure;
plot3(x,y,z,'r-','LineWidth',1);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Full Particle Track in Earth Magnetic Field');
axis equal;
grid on;

%局部回旋轨迹
figure;
plot3(x(1000:6000),y(1000:6000),z(1000:6000),'b-','LineWidth',1);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Zoomed-in Gyration motion');
axis auto;
grid on;
