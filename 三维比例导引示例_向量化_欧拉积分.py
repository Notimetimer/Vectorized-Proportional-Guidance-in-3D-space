import numpy as np
import scipy
import matplotlib.pyplot as plt
from math import pi, sin, cos, tan, atan, atan2

# 有bug, 重力, 推力和阻力对加速度的计算没有使用矢量和

g = 9.81
dt = 0.05
m = 1
kill_range = 5


def Killed(p_m0, p_m1, p_t, kill_range):
    # 直线的两个特征点
    s0 = p_m0  # np.array([1, 0, 0])
    s1 = p_m1  # np.array([-1, 0, 0])
    # 待求解的点
    t0 = p_t  # np.array([0, 0, 3])
    u = s1 - s0
    ct = np.dot(t0 - s0, u) / np.dot(u, u)
    if 0 <= ct <= 1:
        # 点到直线的距离
        dist = np.linalg.norm(s0 + ct * u - t0)
        closest = s0 + ct * u
    else:
        d0 = np.linalg.norm(s0 - t0)
        d1 = np.linalg.norm(s1 - t0)
        if d0 < d1:
            closest = s0
            dist = d0
        else:
            closest = s1
            dist = d1

    if dist <= kill_range:
        return True, closest
    else:
        return False, closest


# 提取最后一行
def latest(vectors):
    if vectors.ndim == 1:  # 一维数组
        return vectors  # 将一维数组转换为二维数组
    elif vectors.ndim == 2:  # 二维数组
        return vectors[-1]  # 直接返回最后一行


# 北天东坐标
pm_ = np.array([-30, 200, 0])
pt_ = np.array([0, 100, 0])

psi_m = 0 * pi / 180  # 导弹初始航向角（自北向东为正）
theta_m = 0 * pi / 180  # 导弹初始俯仰角（自水平面向上为正）
# 导弹不考虑滚转,目标不考虑机动
psi_t = 18 * pi / 180  # 目标初始航向角（自北向东为正）
theta_t = 0 * pi / 180  # 目标初始俯仰角（自水平面向上为正）
vm = 10  # 导弹初始速度大小
vt = 20  # 目标初始速度大小

vm_ = vm * np.array(
    [cos(theta_m) * cos(psi_m),
     sin(theta_m),
     cos(theta_m) * sin(psi_m)])
vt_ = vt * np.array(
    [cos(theta_t) * cos(psi_t),
     sin(theta_t),
     cos(theta_t) * sin(psi_t)])

# Nx = np.array([0])
# Ny = np.array([0])
# Nz = np.array([0])
Nx = []
Ny = []
Nz = []
Vm = []
Line__ = []
Distance = []
killed = False
Q_beta = []
Q_epsilon = []
t_range = np.arange(0, 40, dt)

for t in t_range:
    t_max = t
    pmt_ = latest(pm_)
    ptt_ = latest(pt_)
    vmt_ = latest(vm_)
    vtt_ = latest(vt_)
    line_t_ = ptt_ - pmt_
    distance = np.linalg.norm(line_t_)
    vmt = np.linalg.norm(vmt_)
    vtt = np.linalg.norm(vtt_)

    # 加速参数
    Fx = 0.001 * vmt ** 2  # 阻力
    # Fp = 0
    Fp = 2 * g * m if t < 20 else 0  # 推力

    psi_mt = np.arctan2(vmt_[2], vmt_[0])
    psi_tt = np.arctan2(vtt_[2], vtt_[0])

    vm_hor = np.linalg.norm([vmt_[2], vmt_[0]])  # 导弹水平分速度
    vt_hor = np.linalg.norm([vtt_[2], vtt_[0]])  # 目标水平分速度
    distance_hor = np.linalg.norm([line_t_[2], line_t_[0]])  # 弹目线水平距离

    q_beta_t = np.arctan2(line_t_[2], line_t_[0])  # 目标线偏角
    q_epsilon_t = np.arctan2(line_t_[1], distance_hor)  # 目标线倾角

    theta_mt = np.arctan2(vmt_[1], vm_hor)
    theta_tt = np.arctan2(vtt_[1], vt_hor)

    g_ = g * np.array([0, -1, 0])

    g_parallel_ = np.dot(g_, vmt_) / vmt ** 2 * vmt_

    aT = (Fp - Fx) / m - g * sin(theta_mt)
    aT_ = aT * vmt_ / vmt if vmt > 0 else np.zeros([1, 3])

    # 向量化制导算法
    vrt_ = vtt_ - vmt_
    distance_dot = np.abs(np.dot(vrt_, line_t_) / distance)
    vr_parallel_ = np.dot(vrt_, line_t_) / distance ** 2 * line_t_  # why? 所需过载量垂直于弹目连线
    # vr_parallel_ = np.dot(vrt_, vmt_) / vmt ** 2 * vmt_  # why? 所需过载量垂直于速度
    vr_perpendicular_ = vrt_ - vr_parallel_

    # 参考导引律公式
    # aN_target_required_ = 4 * np.linalg.norm(np.cross(line_t_, vrt_)) / (
    #         distance ** 2) * vmt * vr_perpendicular_ / np.linalg.norm(
    #     vr_perpendicular_) if np.linalg.norm(vr_perpendicular_) > 0 else np.array([0, 0, 0])  # PPG
    # aN_target_required_ = 4 * np.linalg.norm(np.cross(line_t_, vrt_)) / (
    #         distance ** 2) * distance_dot * vr_perpendicular_ / np.linalg.norm(
    #     vr_perpendicular_) if np.linalg.norm(vr_perpendicular_) > 0 else np.array([0, 0, 0])  # TPG

    # # 防止所需侧向过载在速度方向有投影, 处理方式1: 成功但不太合理
    # aN_target_required_ = 4 * np.linalg.norm(np.cross(line_t_, vrt_)) / (
    #         distance ** 2) * distance_dot * vr_perpendicular_ / np.linalg.norm(
    #     vr_perpendicular_) if np.linalg.norm(vr_perpendicular_) > 0 else np.array([0, 0, 0])  # TPG
    # aN_course_required_ = np.zeros([1, 3])  # debug
    # aN_course_required_ = -(g_ - g_parallel_)
    # aN1_ = aN_target_required_ + aN_target_required_
    # aN2 = np.clip(np.linalg.norm(aN1_), 0, 20 * g)
    # aN2_ = aN1_ * aN2 / np.linalg.norm(aN1_) if np.linalg.norm(aN1_) > 0 else np.array([0, 0, 0])
    # aN_ = aN2_ + (g_ - g_parallel_)
    # # test 消去"侧向过载"中平行于速度方向的部分
    # aN_ = aN_ - np.dot(aN_, vmt_) / vmt ** 2 * vmt_

    # # 防止所需侧向过载在速度方向有投影, 处理方式2: 正确
    # aN_target_required_ = 4 * np.linalg.norm(np.cross(line_t_, vrt_)) / (
    #         distance ** 2) * distance_dot * vr_perpendicular_ / np.linalg.norm(
    #     vr_perpendicular_) if np.linalg.norm(vr_perpendicular_) > 0 else np.array([0, 0, 0])  # TPG
    # # aN_course_required_ = np.zeros([1, 3])  # debug
    # aN_course_required_ = -(g_ - g_parallel_)
    # aN1_ = aN_target_required_ + aN_target_required_
    # # test 消去"侧向过载"中平行于速度方向的部分
    # aN1_ = aN1_ - np.dot(aN1_, vmt_) / vmt ** 2 * vmt_
    # aN2 = np.clip(np.linalg.norm(aN1_), 0, 20 * g)
    # aN2_ = aN1_ * aN2 / np.linalg.norm(aN1_) if np.linalg.norm(aN1_) > 0 else np.array([0, 0, 0])
    # aN_ = aN2_ + (g_ - g_parallel_)

    # 防止所需侧向过载在速度方向有投影, 处理方式2.5  正确
    # vrL_ = np.dot(vrt_, line_t_) / distance ** 2 * line_t_
    # vr_perpendicular_ = vrt_ - vrL_
    # n_normalized_ = vr_perpendicular_ / np.linalg.norm(
    #     vr_perpendicular_) if np.linalg.norm(vr_perpendicular_) > 0 else np.array([0, 0, 0])
    # aN_target_required_ = 4 * np.linalg.norm(np.cross(line_t_, vrt_)) / (
    #         distance ** 2) * distance_dot * n_normalized_  # TPG
    # # aN_course_required_ = np.zeros([1, 3])  # debug
    # aN_course_required_ = -(g_ - g_parallel_)
    # aN1_ = aN_target_required_ + aN_target_required_
    # # test 消去"侧向过载"中平行于速度方向的部分
    # aN1_ = aN1_ - np.dot(aN1_, vmt_) / vmt ** 2 * vmt_
    # aN2 = np.clip(np.linalg.norm(aN1_), 0, 20 * g)
    # aN2_ = aN1_ * aN2 / np.linalg.norm(aN1_) if np.linalg.norm(aN1_) > 0 else np.array([0, 0, 0])
    # aN_ = aN2_ + (g_ - g_parallel_)

    # 防止所需侧向过载在速度方向有投影, 处理方式3: 正确
    vrL_ = np.dot(vrt_, line_t_) / distance ** 2 * line_t_  # 径向相对速度
    vr_perpendicular_ = vrt_ - vrL_
    n1_ = vr_perpendicular_  # 周向相对速度
    n2 = n1_ - np.dot(n1_, vmt_) / vmt ** 2 * vmt_  # 周向相对速度在导弹法平面的投影

    n_normalized_ = n2 / np.linalg.norm(n2) if np.linalg.norm(n2) > 0 else np.array([0, 0, 0])

    aN_target_required_ = 4 * np.linalg.norm(np.cross(line_t_, vrt_)) / (
            distance ** 2) * distance_dot * n_normalized_  # TPG
    aN_course_required_ = -(g_ - g_parallel_)
    aN1_ = aN_target_required_ + aN_target_required_
    aN2 = np.clip(np.linalg.norm(aN1_), 0, 20 * g)
    aN2_ = aN1_ * aN2 / np.linalg.norm(aN1_) if np.linalg.norm(aN1_) > 0 else np.array([0, 0, 0])  # 动作引起的加速度
    aN_ = aN2_ + (g_ - g_parallel_)  # 加速度在法向的合量

    # 欧拉积分更新速度
    # vmt += np.dot(aT_ + aN_target_required_ + g * np.array([0, -1, 0]), vmt_) / vmt * dt
    am_ = aT_ + aN_
    vmt += aT * dt

    Transform = np.dot(
        [[cos(theta_mt), sin(theta_mt), 0],
         [-sin(theta_mt), cos(theta_mt), 0],
         [0, 0, 1]],
        [[cos(psi_mt), 0, sin(psi_mt)],
         [0, 1, 0],
         [-sin(psi_mt), 0, cos(psi_mt)]]
    )
    aTNB_ = np.dot(Transform, aN2_)
    Nx.append(aTNB_[0] / g)
    Ny.append(aTNB_[1] / g)
    Nz.append(aTNB_[2] / g)

    # Nx.append(aT / g)
    # Ny.append(np.dot(am_,[0,1,0]) / g)
    # Nz.append(np.dot(am_,[0,0,1]) / g)

    vmt_ += am_ * dt
    vmt_ = vmt_ * vmt / np.linalg.norm(vmt_) if np.linalg.norm(vmt_) > 0 else np.zeros([1, 3])

    if psi_mt > pi:
        psi_mt -= 2 * pi
    if psi_mt < -pi:
        psi_mt += 2 * pi

    # 运动学方程：更新位置
    # vmt_ = vmt * np.array([cos(theta_mt) * cos(psi_mt), sin(theta_mt), cos(theta_mt) * sin(psi_mt)])
    vm_ = np.vstack((vm_, vmt_))
    vt_ = np.vstack((vt_, vtt_))
    pmt_ = pmt_ + vmt_ * dt
    ptt_ = ptt_ + vtt_ * dt
    pm_ = np.vstack((pm_, pmt_))
    pt_ = np.vstack((pt_, ptt_))

    # Ny.append(nyt)
    # Nz.append(nzt)
    Vm.append(vmt)
    Line__.append(line_t_)
    Distance.append(distance)
    Q_beta.append(q_beta_t)
    Q_epsilon.append(q_epsilon_t)

    # 爆炸判断：
    # if distance <= 3:
    #     break
    if t >= 0 + dt:
        killed, point = Killed(pm_[-2], pm_[-1], ptt_, kill_range)
        if killed:
            print('Target killed')
            break

    if vmt <= 0:
        break

# Ny = np.array(Ny).reshape(-1, 1)
# Nz = np.array(Nz).reshape(-1, 1)
Vm = np.array(Vm).reshape(-1, 1)
Line__ = np.array(Line__).reshape(-1, 3)
Distance = np.array(Distance).reshape(-1, 1)
if not killed:
    point = pm_[-1]
# if point == None:
#     point = pm_[-1]

# 画图显示轨迹
p_m_show = np.array(pm_).T
p_t_show = np.array(pt_).T

# 绘制爆炸范围球
u = np.linspace(0, 2 * np.pi, 100)  # 用参数方程画图
v = np.linspace(0, np.pi, 100)
x = point[2] + kill_range * 3 * np.outer(np.cos(u), np.sin(v))  # outer()外积函数：返回cosu × sinv
y = point[0] + kill_range * 3 * np.outer(np.sin(u), np.sin(v))  #
z = point[1] + kill_range * 3 * np.outer(np.ones(np.size(u)), np.cos(v))  # ones()：返回一组[1,1,.......]

# plt.ion()
fig = plt.figure(1)
ax3d = fig.add_subplot(projection='3d')
ax3d.plot(p_m_show[2], p_m_show[0], p_m_show[1], c='b', label='MissleTrack')  # 显示为东北天
ax3d.scatter(p_m_show[2][0], p_m_show[0][0], p_m_show[1][0], c='b', s=5)
ax3d.plot(p_t_show[2], p_t_show[0], p_t_show[1], c='r', label='TargetTrack')
ax3d.scatter(p_t_show[2][0], p_t_show[0][0], p_t_show[1][0], c='r', s=5)
# ax3d.scatter(point[2], point[0], point[1], c='pink', s=kill_range)
ax3d.plot_surface(x, y, z, cmap='hot', alpha=0.3)  # hot 或是 cool，其他没试过


# 设置坐标轴等比例
def set_axes_equal(ax):
    """确保3D图的坐标轴单位长度相等。"""
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max(x_range, y_range, z_range)
    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)
    ax.set_xlim3d([x_middle - max_range / 2, x_middle + max_range / 2])
    ax.set_ylim3d([y_middle - max_range / 2, y_middle + max_range / 2])
    ax.set_zlim3d([z_middle - max_range / 2, z_middle + max_range / 2])
    # 东北天
    ax.set_xlabel('Z')
    ax.set_ylabel('X')
    ax.set_zlabel('Y')
    # 设置等显示缩放比例
    ax.set_box_aspect([1, 1, 1])


# 调用函数设置等比例坐标轴
set_axes_equal(ax3d)
# ax3d.set_xlim3d(left=-1000, right=1000)  # 设置x轴的范围从0到6
# ax3d.set_ylim3d(bottom=-1000, top=1000)  # 设置y轴的范围从0到6
# ax3d.set_zlim3d(bottom=-1000, top=1000)  # 设置z轴的范围从0到6

'''
补充弹目距离, 弹目线俯仰角, 弹目线偏角, 纵向与侧向过载量, 导弹速度关于时间的曲线
'''
plt.figure(2)
plt.subplot(411)
plt.plot(range(len(Distance)), Distance, color='r', linestyle='-', marker='', label='distance')
plt.legend(loc='upper right')
plt.subplot(412)
plt.plot(range(len(Distance)), Q_beta, color='r', linestyle='-', marker='', label='q_beta')
plt.plot(range(len(Distance)), Q_epsilon, color='b', linestyle='-', marker='', label='q_epsilon')
plt.legend(loc='upper right')
plt.subplot(413)
plt.plot(range(len(Distance)), Nx, color='r', linestyle='-', marker='', label='nx')
plt.plot(range(len(Distance)), Ny, color='g', linestyle='-', marker='', label='ny')
plt.plot(range(len(Distance)), Nz, color='b', linestyle='-', marker='', label='nz')
plt.legend(loc='upper right')
plt.subplot(414)
plt.plot(range(len(Distance)), Vm, color='r', linestyle='-', marker='', label='vm')
plt.legend(loc='upper right')
plt.show()
# plt.show(block=True)
