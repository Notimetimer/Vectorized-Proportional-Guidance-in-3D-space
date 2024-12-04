clc
clear
close all

dist0=5e3;
heightt0=10e3;
psi_m=200*pi/180;
theta_m=-90*pi/180;
theta_t=0;
vm=360;
vt=360;
shoot_method=0;
target_move=0;
count=0;
% [hit,count]=MissileTest(dist0, heightt0, psi_m, theta_m, theta_t, vm, vt,shoot_method, target_move, count)
% function [hit,count]=MissileTest(dist0, heightt0, psi_m, theta_m, theta_t, vm, vt,shoot_method, target_move, count)
draw_curve=1;
hit=0;
% dist0=3e3;
% heightt0=5000+-200;
% psi_m = 160 * pi / 180; % initial missile heading angle (positive from north to east)
% theta_m = -80 * pi / 180; % initial missile pitch angle (positive from horizontal to up)
% vm = 150; % initial missile speed
% vt = 100; % initial target speed
% theta_t = 0 * pi / 180; % initial target pitch angle (positive from horizontal to up)

% Initialize constants
g = 9.81;
dt = 0.02;
m = 87; % kg
length = 2.87; % m (missile length)
Diameter = 0.127; % m (missile diameter)
S = pi * (Diameter / 2)^2; % missile cross-sectional area
dm = 6; % fuel flow rate per second
Isp = 170; % specific impulse
kill_range = 15;

psi_t = 0 * pi / 180; % initial target heading angle (positive from north to east)

% Initial positions and velocities
pt_ = [0, heightt0, 0];

heightm0=heightt0-dist0*sin(theta_m);

pm_ = [pt_(1)-dist0*cos(theta_m)*cos(psi_m), ...
    heightm0, ...
    pt_(3)-dist0*cos(theta_m)*sin(psi_m)];


% 平射/对准目标发射
if shoot_method==0
    theta_m=min(theta_m, pi/2);
elseif shoot_method==1
    theta_m=0; 
end
theta_m=theta_m(1);

vm_ = vm * [cos(theta_m) * cos(psi_m), sin(theta_m), cos(theta_m) * sin(psi_m)];
vt_ = vt * [cos(theta_t) * cos(psi_t), sin(theta_t), cos(theta_t) * sin(psi_t)];

% vm_ = vm * [cos(theta_m) * cos(psi_m); sin(theta_m); cos(theta_m) * sin(psi_m)];
% vt_ = vt * [cos(theta_t) * cos(psi_t); sin(theta_t); cos(theta_t) * sin(psi_t)];

Nx = [];
Ny = [];
Nz = [];
Vm = [];
Line__ = [];
Distance = [];
killed = false;
Q_beta = [];
Q_epsilon = [];
t_range = 0:dt:60;

last_pmt_=pm_;
last_ptt_=pt_;
last_vmt_=vm_;
last_vtt_=vt_;

for t = t_range
%     t_max = t;
    if draw_curve==1
        pmt_ = latest(pm_);
        ptt_ = latest(pt_);
        vmt_ = latest(vm_);
        vtt_ = latest(vt_);
    else
        pmt_=last_pmt_;
        ptt_=last_ptt_;
        vmt_=last_vmt_;
        vtt_=last_vtt_;
    end

    line_t_ = ptt_ - pmt_;
    distance = norm(line_t_);
    vmt = norm(vmt_);
    vtt = norm(vtt_);

    % 导弹动力学、运动学与导引率
    % 导弹马赫数
    height=pmt_(2);
    if height<=11000
        soundspeed=20.05*sqrt(288.15-0.00651122*height);
    else
        soundspeed=295.069;
    end
    ma=vmt/soundspeed;
    % 导弹阻力系数
    if ma<0.8
        Cd=0.5;
    elseif ma>=0.8 && ma<1.2
        Cd=0.5+(ma-0.8)*(0.78-0.5)/0.4;
    elseif ma>=1.2 && ma<3
        Cd=1.23582-0.495742*ma+0.108388*ma^2+-0.00979703*ma^3;
    else 
        Cd=1;
    end

    % Acceleration parameters
    rho = 1.225 * exp(-pmt_(2) / 9300); % air density in kg/m^3
    Fx = 0.5 * Cd * S * rho * vmt^2; % drag force
    Fp = g * Isp * dm * (t < 5.2); % thrust force

    psi_mt = atan2(vmt_(3), vmt_(1));
    if t<dt
        psi_tt = atan2(vtt_(3), vtt_(1)); % debug
    end
    vm_hor = norm([vmt_(3), vmt_(1)]); % missile horizontal velocity
    vt_hor = norm([vtt_(3), vtt_(1)]); % target horizontal velocity
    distance_hor = norm([line_t_(3), line_t_(1)]); % horizontal distance between missile and target

    q_beta_t = atan2(line_t_(3), line_t_(1)); % target line azimuth angle
    q_epsilon_t = atan2(line_t_(2), distance_hor); % target line elevation angle

    theta_mt = atan2(vmt_(2), vm_hor);
    theta_tt = atan2(vtt_(2), vt_hor);

    g_ = g * [0, -1, 0];
    g_parallel_ = dot(g_, vmt_) / vmt^2 * vmt_;

    if t <5.2
        m = m - dm * dt;
    end
    aT = (Fp - Fx) / m - g * sin(theta_mt);
    aT_ = aT * vmt_ / vmt;

    % Vectorized guidance algorithm
    vrt_ = vtt_ - vmt_;
    distance_dot = abs(dot(vrt_, line_t_) / distance);
    vr_parallel_ = dot(vrt_, line_t_) / distance^2 * line_t_;
%     vr_perpendicular_ = vrt_ - vr_parallel_;

    vrL_ = dot(vrt_, line_t_) / distance^2 * line_t_;
    vr_perpendicular_ = vrt_ - vrL_;
    n1_ = vr_perpendicular_;
    n2_ = n1_ - dot(n1_, vmt_) / vmt^2 * vmt_;
    if norm(n2_)>0
        n_normalized_ = n2_ / norm(n2_);
    else
        n_normalized_=[0,0,0];
    end

    aN_target_required_ = 4 * norm(cross(line_t_, vrt_)) / (distance^2) * distance_dot * n_normalized_;
        
    % 导引头角度限制
    if dot(vmt_,line_t_)/(vmt*distance)<cos(90*pi/180) % 范围+-90°
        aN_target_required_=[0,0,0];
    end
    % 导引头距离限制
    if distance>10e3
        aN_target_required_=[0,0,0];
    end
    % 导引头角速度限制
    if norm(vr_perpendicular_)/distance>28*pi/180
        aN_target_required_=[0,0,0];
    end

    aN_course_required_ = -(g_ - g_parallel_);

    aN_course_required_ = aN_course_required_...
        -dot(aN_course_required_, vmt_)/ vmt^2 * vmt_;
    
    aN1_ = aN_target_required_ + aN_course_required_;
    aN2 = min(max(norm(aN1_), 0), 30 * g);
    if norm(aN1_)>0
        aN2_ = aN1_ * aN2 / norm(aN1_);
    else
        aN2_ = [0,0,0];
    end
    aN_ = aN2_ + (g_ - g_parallel_);

    % Euler integration to update velocity
    am_ = aT_ + aN_;
    vmt = vmt + aT * dt;

    Transform = [cos(theta_mt), sin(theta_mt), 0; -sin(theta_mt), cos(theta_mt), 0; 0, 0, 1] * ...
                [cos(psi_mt), 0, sin(psi_mt); 0, 1, 0; -sin(psi_mt), 0, cos(psi_mt)];
    aTNB_ = (Transform * aN2_')';

    if draw_curve==1
        Nx = [Nx, aTNB_(1) / g];
        Ny = [Ny, aTNB_(2) / g];
        Nz = [Nz, aTNB_(3) / g];
    end
    vmt_ = vmt_ + am_ * dt;
    % 限速
    if vmt/soundspeed>3
        vmt=soundspeed*3;
    end
    vmt_ = vmt_ * vmt / norm(vmt_);
    
    psi_mt=limit_angle(psi_mt);

    pmt_ = pmt_ + vmt_ * dt;

    %%
    % 目标运动学
    % todo 目标逃逸模型放在这里
    % Simulation parameters
%     Missile = HighProjectileTrajectory();  % 实例化类
%     changelog：目标运动状态更新
    % 目标运动状态更新
    
    % 在此规定逃逸角分为水平和垂直，是由目标线转向目标速度矢量的角度，水平以又转为正，垂直以向上为正
    % 水平逃逸角写作EAh，垂直逃逸角度写为EAv

    % 垂直加速度参考值
    ay_refer = 6*g;
    % 水平加速度参考值
    az_refer = 6*g;
    EAh=psi_tt-q_beta_t;
    EAv=theta_tt-q_epsilon_t;
    if EAh > pi
        EAh = EAh - 2 * pi;
    end
    if EAh <= -pi
        EAh = EAh + 2 * pi;
    end

    % 0 平飞
    if target_move==0
        ay_t=0;
        az_t=0;
    end
    % 1 39线机动
    if target_move==1
        ay_t=0;
        if 0<=EAh&&EAh<pi/2 || -pi<EAh&&EAh<-pi/2 % 如果被导弹尾追或是拦堵
            az_t=az_refer;
        elseif -pi/2<EAh&&EAh<0 || pi/2<EAh&&EAh<=pi
            az_t=-az_refer;
        else
            az_t=0;
        end
%         disp(az_t) % debug
    end
    % 2 水平置尾机动
    if target_move==2
        ay_t=0;
%         disp(EAh*180/pi)
        if 0<EAh&&EAh<=pi % 0<EAh&&EAh<=pi
            az_t=-az_refer;
        elseif -pi<EAh&&EAh<0
            az_t=az_refer;
        else
            az_t=0;
        end
    end
    %{
        目标机动起始设定：当距离略大于导弹末制导距离时目标才可以开始做机动
        另一种机动起始点设定方法：计算角速度和导弹剩余命中时间，并保证剩余命中时间刚够转pi/2
        爬升和俯冲机动的起始时刻必须
    %}
    % 3 爬升机动
    if target_move==3
        az_t=0;
        if distance<=10e3*2 %Missile.detect_distance*1
            ay_t=ay_refer;
        else
            ay_t=0;
        end
%         if ptt_(2)<15000
%             ay_t=ay_refer;
%         else
%             ay_t=ay_refer/exp((ptt_(2)-15000)/600)*-sin(theta_tt);
%         end
    end
    % 4 俯冲机动
    if target_move==4
        az_t=0;
        if distance<=10e3*2 %Missile.detect_distance*1
            ay_t=-ay_refer;
        else
            ay_t=0;
        end
%         if ptt_(2)>3000
%             ay_t=-ay_refer;
%         else
%             ay_t=0;
%             theta_tt=0;
% %             ay_t=ay_refer/exp((ptt_(2)-2000)/200)*-sin(theta_tt);
%         end
    end

%     % 5 水平蛇形机动
%     if target_move==5
%         ay_t=0;
%         az_t=az_refer*cos(2*pi* 0.06 *t);
%     end
    % changelog

    % 5 水平蛇形机动(形状相似)
    if target_move==5
        ay_t=0;
        az_t=az_refer*cos(2*pi* 0.06 * t * 0.5*340/vtt);
    end

    % 动力学模型更新

    % Euler angle anti-singularity
    theta_limit=85*pi/180; %Missile.theta_limit;

    if height <= 11000
    sound_speed = 20.0463 * sqrt(288.15 - 0.00651122 * height);
    else
    sound_speed = 295.069;
    end
       
    % 动力学模型更新
    vtt=vtt+0*2*g*dt;
    vtt=min(vtt,1.8*sound_speed); % debug
    theta_tt = theta_tt + (ay_t/ vtt) * dt;
% debug
%     if ay_t<0
%         theta_tt=-30*pi/180;
%     end
    psi_tt = psi_tt + az_t/vtt/cos(theta_tt)*dt;
    theta_tt = max(min(theta_tt, theta_limit), -theta_limit); % debug
    theta_tt=theta_tt(1);
    % debug
    if ptt_(2)<3000 || ptt_(2)>18000
        theta_tt=0;
    end

    if psi_tt > pi
        psi_tt = psi_tt - 2 * pi;
    end
    if psi_tt <= -pi
        psi_tt = psi_tt + 2 * pi;
    end
    vtt_ = vtt * [cos(theta_tt) * cos(psi_tt), sin(theta_tt), cos(theta_tt) * sin(psi_tt)];
%   欧拉积分更新位置
    ptt_ = ptt_ + vtt_ * dt;

%%
    % 伤害判定
    if t >= 0 + dt
        [killed, point] = Killed(0, pmt_, ptt_, kill_range);
        if killed && norm(vrt_)>=30 % 最小引信速度
%             disp('Target killed');
            hit=1;
            break;
        end
    end

    if vmt <= 50
        break;
    end
    if pmt_(2)<100
        break;
    end
    if t>5.2*2
        if distance>last_distance
            break;
        end
    end

    if draw_curve==1
        % Kinematic equations: update positions
        vm_ = [vm_; vmt_];
        vt_ = [vt_; vtt_];
        pm_ = [pm_; pmt_];
        pt_ = [pt_; ptt_];
        Vm = [Vm, vmt];
        Line__ = [Line__; line_t_];
        Distance = [Distance; distance];
        Q_beta = [Q_beta; q_beta_t];
        Q_epsilon = [Q_epsilon; q_epsilon_t];
        last_distance=distance;
    else
        last_pmt_=pmt_;
        last_ptt_=ptt_;
        last_vmt_=vmt_;
        last_vtt_=vtt_;
        last_distance=distance;
    end

end

Vm = Vm';
Line__ = Line__';
Distance = Distance';
if ~killed
    point = pm_(end, :);
end


if draw_curve==1
    % Plot trajectories
    p_m_show = pm_';
    p_t_show = pt_';

    figure(1);
    hold off;
    plot3(p_m_show(3, :), p_m_show(1, :), p_m_show(2, :), 'b', 'DisplayName', 'MissileTrack');
    hold on;
    scatter3(p_m_show(3, 1), p_m_show(1, 1), p_m_show(2, 1), 'b', 'filled');
    plot3(p_t_show(3, :), p_t_show(1, :), p_t_show(2, :), 'r', 'DisplayName', 'TargetTrack');
    scatter3(p_t_show(3, 1), p_t_show(1, 1), p_t_show(2, 1), 'r', 'filled');
    % sphere画爆炸火球

    % Set equal axis scaling
    axis equal;
    xlabel('Z');
    ylabel('X');
    zlabel('Y');

    % Additional plots
    figure(2);
    subplot(4, 1, 1);
    plot(Distance, 'r');
    title('Distance');
    subplot(4, 1, 2);
    plot(Q_beta, 'r');
    hold on;
    plot(Q_epsilon, 'b');
    hold off;
    title('Q_{beta} and Q_{epsilon}');
    legend('Q_{beta}','Q_{epsilon}')
    subplot(4, 1, 3);
    plot(Nx, 'r');
    hold on;
    plot(Ny, 'g');
    plot(Nz, 'b');
    hold off;
    title('Nx, Ny, and Nz');
    legend('Nx','Ny','Nz')
    subplot(4, 1, 4);
    plot(Vm, 'r');
    title('Vm');
end

% Function to check if the target is killed
function [killed, closest] = Killed(~, p_m1, p_t, kill_range)
%     % 求点和线段的距离
%     s0 = p_m0;
%     s1 = p_m1;
%     t0 = p_t;
%     u = s1 - s0;
%     ct = dot(t0 - s0, u) / dot(u, u);
%     if 0 <= ct && ct <= 1
%         dist = norm(s0 + ct * u - t0);
%         closest = s0 + ct * u;
%     else
%         d0 = norm(s0 - t0);
%         d1 = norm(s1 - t0);
%         if d0 < d1
%             closest = s0;
%             dist = d0;
%         else
%             closest = s1;
%             dist = d1;
%         end
%     end

%     % 求点和点距离
    dist=norm(p_m1-p_t);
    closest=p_m1;

    if dist <= kill_range
        killed = true;
    else
        killed = false;
    end
end

% Function to get the last row of a matrix
function last_row = latest(vectors)
    if size(vectors, 1) == 1
        last_row = vectors;
    else
        last_row = vectors(end, :);
    end
end

function out=limit_angle(in)
    if in<=-pi
        out=in+2*pi;
    elseif in>pi
        out=in-2*pi;
    else
        out=in;
    end
end