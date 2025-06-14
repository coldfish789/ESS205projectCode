close all;clear all;clc

fid = fopen("track_RK.bin",'rb');
data = fread(fid,[3,Inf],'double');
fclose(fid);

x = data(1,:);
y = data(2,:);
z = data(3,:);

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


% 投影到 XY 平面
figure('Position', [100, 100, 1200, 400]);
subplot(1,3,1);
plot(x, y, "r-", "LineWidth", 1);
xlabel('X');
ylabel('Y');
title('XY Projection (Drift Motion)');
axis equal;
grid on;

% 投影到 XZ 平面
subplot(1,3,2);
plot(x, z, "g-", "LineWidth", 1);
xlabel('X');
ylabel('Z');
title('XZ Projection(bounce)');
grid on;

% 投影到 YZ 平面
subplot(1,3,3);
plot(y, z, "b-", "LineWidth", 1);
xlabel('Y');
ylabel('Z');
title('YZ Projection(bounce)');
grid on;
