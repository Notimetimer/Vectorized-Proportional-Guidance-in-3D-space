clear
close all
clc

d_epsilon= 10 * pi/180; %20
d_beta= 20 * pi/180;
R_test=50;
count=0;

Q_epsilon=-pi/2:d_epsilon:pi/2;
Q_beta=-pi+d_beta:d_beta:pi;

max_height=20e3;
min_height=1e3;

% 搜索超参数
a0=800;
b0=30e3;
b0_close=10e3;
min_range=400;

epsilon = 500; % 搜索精度
delta = 5000;

Dot_=[];
heightt0=10000;
pt0_=[0,heightt0,0]; % 目标初始位置
vm0 = 560; % initial missile speed
vt0 = 560; % initial target speed
theta_t = 0 * pi / 180; % initial target pitch angle (positive from horizontal to up)
% 发射方式
shoot_method=0;
%{
0：瞄准目标发射
1：平射
2：仰射（待做）
%}
% 目标机动样式
target_move=0;
%{
目标机动编号
0、平飞
1、39线机动
2、水平置尾机动
3、爬升机动
4、俯冲机动
5、水平蛇形机动
6、1g过载右转
%}
% dist0=b0
% q_epsilon=Q_epsilon(1)
% q_beta=Q_beta(1)
% [hit,~]=AIM9Ltest_2(dist0, heightt0, q_beta, q_epsilon, theta_t, vm0, vt0, shoot_method, target_move, count);

angle_list=[];
for q_epsilon=Q_epsilon % -Q_epsilon debug
    for q_beta=Q_beta
        % 筛选掉南北极多余的经度
        if abs(q_epsilon)==pi/2 && q_beta~=Q_beta(1)
            continue
        end
        % 临时存储仿真结果
        angle_list=[angle_list;[q_epsilon,q_beta]];
    end
end
L=size(angle_list,1);

run_list=zeros(L,2);
% note 跳过蒙特卡洛直接搜索边界，如果发现错误需要回到蒙特卡洛去

% i=18
parfor i = 1:L % parfor
    q_epsilon=angle_list(i,1);
    q_beta=angle_list(i,2);
    % ↓在这里计算某个初始状态下的攻击远近边界，放入run_list中
    % 计算远边界
    a = a0;
    if q_epsilon<0
        temp=(max_height-heightt0)/(abs(sin(q_epsilon))+1e-10);
        b=min([b0, temp]);
        b=b(1);
    elseif q_epsilon>0
        temp=(heightt0-min_height)/(abs(sin(q_epsilon))+1e-10);
        b = min([b0, temp]);
        b=b(1);
    else
        temp=100e3; % 设置一个不可能达到的大值
        b=b0;
    end
    if b<b0
        go_farther=0; % 不往高度限制外搜索
    else
        go_farther=1;
    end
    next_hit1=-1;
    next_hit2=-1;
    next_hit3=-1;
    iter=0;
    count=0;
    while abs(a-b)>epsilon
        iter=iter+1;
        lambda=a+0.5*(b-a);
        % 远边界不用算a点
        % 算λ点
        dist0=lambda;
        [hit2,count]=AIM9Ltest_2(dist0, heightt0, q_beta, q_epsilon, theta_t, vm0, vt0, shoot_method, target_move, count);
        % 算b点
        if next_hit3==-1
            dist0=b;
            [hit3,count]=AIM9Ltest_2(dist0, heightt0, q_beta, q_epsilon, theta_t, vm0, vt0, shoot_method, target_move, count);
        else
            hit3=next_hit3;
        end
        next_hit3=-1;

        if hit3==1 && go_farther==1
            b=b+delta;
            continue
        end
        if hit2==1
            a=lambda;
            next_hit1=hit2;
        else
            b=lambda;
            next_hit3=hit2;
        end
    end
    far_range = (a+b)/2;
    % 远边界限高
    if q_epsilon<0
        temp=(max_height-heightt0)/(abs(sin(q_epsilon))+1e-10);
        far_range=min([far_range, temp]);
        far_range=far_range(1);
    elseif q_epsilon>0
        temp=(heightt0-min_height)/(abs(sin(q_epsilon))+1e-10);
        far_range = min([far_range, temp]);
        b=b(1);
    else
        % do nothing
    end

    % 计算近边界
    a = a0;
    if q_epsilon>0
        temp=(max_height-heightt0)/(abs(sin(q_epsilon))+1e-10);
        b=min([b0_close, temp]);
        b = min(b,far_range);
        b=b(1);
    elseif q_epsilon<0
        temp=(heightt0-min_height)/(abs(sin(q_epsilon))+1e-10);
        b = min([b0_close,temp]);
        b = min(b,far_range);
        b=b(1);
    else
        b = min(b0_close,far_range);
        b=b(1);
    end
%     disp('far done')
    % 近边界搜索不会跑出高度限制
%     if b<b0_close
%         go_farther=0; % 不往高度限制外搜索
%     else
%         go_farther=1;
%     end
    count=0;
    next_hit1=-1;
    next_hit2=-1;
    next_hit3=-1;
    iter=0;

    while abs(a-b)>epsilon/5
        iter=iter+1;
        % debug
%         if iter>300
%             break
%         end

        lambda=a+0.5*(b-a);
        % 算a点
        if next_hit1==-1
            dist0=a;
            [hit1,count]=AIM9Ltest_2(dist0, heightt0, q_beta, q_epsilon, theta_t, vm0, vt0, shoot_method, target_move, count);
        else
            hit1=next_hit1;
        end
        % 算λ点
        dist0=lambda;
        [hit2,count]=AIM9Ltest_2(dist0, heightt0, q_beta, q_epsilon, theta_t, vm0, vt0, shoot_method, target_move, count);
        % 近边界搜索不用算b点

        next_hit1=-1;
        next_hit3=-1;

        if hit1==1
            % debug
            if a<=min_range
                b=a; % 0距离能够命中，但是不允许这么打
                break
            end
            a=max(min_range, a-delta/10);
            a=a(1);
            continue;
        end

        if hit2==1
            b=lambda;
            next_hit3=hit2;
        else
            a=lambda;
            next_hit1=hit2;
        end
    end
%     disp('close done')
    close_range = (a+b)/2;

    % 临时存储仿真结果
    run_list(i,:)=[far_range, close_range];
end

%% 
edge_list=[angle_list,run_list];
% edge_list(1,1)=-0.4; % debug
%%
q_epsilon_list=edge_list(:,1);
q_beta_list=edge_list(:,2);
r_far_list=edge_list(:,3);
r_close_list=edge_list(:,4);
% Dot_far=pt0_+r_far_list.*[cos(-q_epsilon_list).*cos(-q_beta_list),sin(-q_epsilon_list),cos(-q_epsilon_list).*sin(-q_beta_list)]; % 仍然使用北天东
Dot_far=pt0_-r_far_list.*[cos(q_epsilon_list).*cos(q_beta_list),sin(q_epsilon_list),cos(q_epsilon_list).*sin(q_beta_list)]; % 仍然使用北天东
% Dot_close=pt0_+r_close_list.*[cos(-q_epsilon_list).*cos(-q_beta_list),sin(-q_epsilon_list),cos(-q_epsilon_list).*sin(-q_beta_list)]; % 仍然使用北天东
Dot_close=pt0_-r_close_list.*[cos(q_epsilon_list).*cos(q_beta_list),sin(q_epsilon_list),cos(q_epsilon_list).*sin(q_beta_list)]; % 仍然使用北天东
far_range_show=show_loc(Dot_far);
close_range_show=show_loc(Dot_close);

% figure(1)
% hold off
% scatter3(edge_list(:,1),edge_list(:,2),edge_list(:,3),'filled','b')
% hold on
% scatter3(edge_list(:,1),edge_list(:,2),edge_list(:,4),'filled','r')
% xlabel('q_{epsilon}');
% ylabel('q_{beta}');
% zlabel('distance');

figure(2)
hold off
% 画远边界点
% scatter3(far_range_show(:,1),far_range_show(:,2),far_range_show(:,3),'filled','b')
hold on
% 画近边界点
% scatter3(close_range_show(:,1),close_range_show(:,2),close_range_show(:,3),'filled','r')
% 围成曲面
far_edge=[edge_list(:,1:3)];
close_edge=[edge_list(:,1:2),edge_list(:,4)];
% 输入必须使用theta慢循环,psi快循环的形式
far_delta=spherical_zone_visualize(pt0_,far_edge,d_epsilon,d_beta,[0.1, 0.3, 0.7],0.15);
close_delta=spherical_zone_visualize(pt0_,close_edge,d_epsilon,d_beta,[0.9, 0.6, 0.1],0.5);
view([-45, 30]);
grid on
% 绘制目标初速指示箭头
% 远边界 箭头
arrow_far_length=max(edge_list(:,3))/2;
arrow_far_length=arrow_far_length(1);
startPoint = [0, 0, heightt0]; % 起点坐标
endPoint = [0, arrow_far_length, heightt0];   % 终点坐标
direction = endPoint - startPoint; % 计算箭头的方向向量
quiver3(startPoint(1), startPoint(2), startPoint(3), ...
         direction(1), direction(2), direction(3), ...
         0, 'LineWidth', 2, 'MaxHeadSize', 1, 'Color', 'r'); % 0表示箭头不缩放
axis("equal")

figure(3)
hold off
close_delta=spherical_zone_visualize(pt0_,close_edge,d_epsilon,d_beta,[0.9, 0.6, 0.1],0.4);
hold on
view([-45, 30]);
grid on
% 近边界箭头
arrow_close_length=max(edge_list(:,4))*1.5;
arrow_close_length=arrow_close_length(1);
endPoint = [0, arrow_close_length, heightt0];
direction = endPoint - startPoint;
quiver3(startPoint(1), startPoint(2), startPoint(3), ...
         direction(1), direction(2), direction(3), ...
         0, 'LineWidth', 2, 'MaxHeadSize', 1, 'Color', 'r'); % 0表示箭头不缩放
scatter3(startPoint(1), startPoint(2), startPoint(3),'filled','b')
axis("equal")


%%
% 东北天转自北天东
function P_show=show_loc(P_)
    P_show=zeros(size(P_));
    P_show(:,1)=P_(:,3);
    P_show(:,2)=P_(:,1);
    P_show(:,3)=P_(:,2);
end