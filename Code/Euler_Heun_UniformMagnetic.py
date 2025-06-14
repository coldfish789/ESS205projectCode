import numpy as np
import matplotlib.pyplot as plt

# 粒子参数
q = 1.0      # 电荷量 (C)
m = 1.0      # 质量 (kg)
Bz = 1.0     # 磁场强度 (T)
B = np.array([0, 0, Bz])

# 初始条件
r0 = np.array([0.0, 0.0, 0.0])
v0 = np.array([1.0, 1.0, 2.0])

# 时间参数
dt = 0.1
T = 3000  # 总时间 (s)
steps = int(T / dt)

# 角频率
omega = q * Bz / m

# 初始化数组
r_euler = np.zeros((steps, 3))
v_euler = np.zeros((steps, 3))
r_cromer = np.zeros((steps, 3))
v_cromer = np.zeros((steps, 3))
r_heun = np.zeros((steps, 3))
v_heun = np.zeros((steps, 3))
r_analytic = np.zeros((steps, 3))

r_euler[0] = r0
v_euler[0] = v0
r_cromer[0] = r0
v_cromer[0] = v0
r_heun[0] = r0
v_heun[0] = v0

# 欧拉法
for i in range(steps-1):
    F = q * np.cross(v_euler[i], B)
    a = F / m
    v_euler[i+1] = v_euler[i] + a * dt
    r_euler[i+1] = r_euler[i] + v_euler[i] * dt

# 梯形欧拉法
for i in range(steps-1):
    # 预测
    F1 = q * np.cross(v_heun[i], B)
    a1 = F1 / m
    v_predict = v_heun[i] + a1 * dt
    r_predict = r_heun[i] + v_heun[i] * dt

    F2 = q * np.cross(v_predict, B)
    a2 = F2 / m

    v_heun[i+1] = v_heun[i] + 0.5 * dt * (a1 + a2)
    r_heun[i+1] = r_heun[i] + 0.5 * dt * (v_heun[i] + v_heun[i+1])

# Analytical solution
t_arr = np.linspace(0, T, steps)
vx0, vy0, vz0 = v0
x_analytic = (vy0/omega) - (vy0/omega)*np.cos(omega*t_arr) + (vx0/omega)*np.sin(omega*t_arr)
y_analytic = -(vx0/omega) + (vx0/omega)*np.cos(omega*t_arr) + (vy0/omega)*np.sin(omega*t_arr)
z_analytic = vz0 * t_arr
r_analytic[:,0] = x_analytic
r_analytic[:,1] = y_analytic
r_analytic[:,2] = z_analytic

# 误差计算：末端位置的欧几里得距离
err_euler = np.linalg.norm(r_euler[-1] - r_analytic[-1])
err_heun = np.linalg.norm(r_heun[-1] - r_analytic[-1])

print(f"Euler method: {err_euler:.6f}")
print(f"Improved Euler method: {err_heun:.6f}")

# 可视化
fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(111, projection='3d')
ax.plot(r_analytic[:,0], r_analytic[:,1], r_analytic[:,2], label='Analytical solution', lw=3)
ax.plot(r_euler[:,0], r_euler[:,1], r_euler[:,2], '--', label='Euler method')
ax.plot(r_heun[:,0], r_heun[:,1], r_heun[:,2], '-.', label='Improved Euler method')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Analytical solution vs Numerical Methods')
ax.legend()
plt.show()

# 误差数组初始化
err_euler_arr = np.zeros(steps)
err_heun_arr = np.zeros(steps)

# 每一步计算误差
for i in range(steps):
    err_euler_arr[i] = np.linalg.norm(r_euler[i] - r_analytic[i])
    err_heun_arr[i] = np.linalg.norm(r_heun[i] - r_analytic[i])

# 绘制误差随时间变化的曲线
plt.figure(figsize=(8,5))
plt.plot(t_arr, err_heun_arr, label="Heun method")
plt.xlabel("Time (s)")
plt.ylabel("Error (Euclidean distance)")
plt.title("Numerical Error vs Analytical Solution")
plt.legend()
plt.grid()
plt.show()

# 计算每一步的 r 和 v 误差
err_r_heun = np.linalg.norm(r_heun - r_analytic, axis=1)
err_v_heun = np.linalg.norm(v_heun - v0, axis=1)  # 解析解速度模长始终与v0相同（磁场不变）

# 可视化
plt.figure(figsize=(10, 5))
plt.plot(t_arr, err_r_heun, label="Position error |r_heun - r_analytic|")
plt.plot(t_arr, err_v_heun, label="Velocity error |v_heun - v_analytic|")
plt.xlabel("Time (s)")
plt.ylabel("Error (Euclidean distance)")
plt.title("Heun Method: Position and Velocity Error vs Analytical Solution")
plt.legend()
plt.grid()
plt.show()

