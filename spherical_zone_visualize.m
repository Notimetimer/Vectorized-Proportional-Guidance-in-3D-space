function [points_delta_xyz]=spherical_zone_visualize(pt0_,one_edge,d_epsilon,d_beta,RGB,transparency)
% 可视化必须容忍重复
% 参数方程造个球
% 创建极角和方位角的采样点
psi = -pi+d_beta:d_beta:pi; % linspace(0, 2*pi, n+1); % 方位角 [0, 2π]
theta = -pi/2:d_epsilon:pi/2; % [-pi/2, π/2]
psi=[psi,-pi+d_beta];

% 创建网格
[thetaGrid, psiGrid] = meshgrid(theta, psi);
alpha=thetaGrid(:);
beta=psiGrid(:);
r_temp=ones(size(psiGrid));
rho=r_temp(:);
points_sph=[alpha,beta,rho];

m=length(theta); % 7
n=length(psi); % 13

% 小的矩阵扩展为中的
temp=[ones(n-2,1)*one_edge(1,:);one_edge;ones(n-2,1)*one_edge(end,:)];
% 中的矩阵扩展为大的
for i=1:m
    points_sph((i-1)*n+(1:n),:)=[temp((i-1)*(n-1)+(1:n-1),:);temp((i-1)*(n-1)+1,:)];
end

alpha_reshaped=reshape(points_sph(:,1),size(thetaGrid));
beta_reshaped=reshape(points_sph(:,2),size(psiGrid));
rho_reshaped=reshape(points_sph(:,3),size(r_temp));

% 球面方程
% x = pt0_(1)+rho_reshaped.*cos(-alpha_reshaped) .* cos(-beta_reshaped);
% z = pt0_(3)+rho_reshaped.*cos(-alpha_reshaped) .* sin(-beta_reshaped);
% y = pt0_(2)+rho_reshaped.*sin(-alpha_reshaped);
x = pt0_(1)-rho_reshaped.*cos(alpha_reshaped) .* cos(beta_reshaped);
z = pt0_(3)-rho_reshaped.*cos(alpha_reshaped) .* sin(beta_reshaped);
y = pt0_(2)-rho_reshaped.*sin(alpha_reshaped);
% points_delta_xyz={x,y,z}; % 北天东格式
points_delta_xyz={z,x,y}; % 东北天格式
% figure()
h=surf(points_delta_xyz{1}, points_delta_xyz{2}, points_delta_xyz{3}); % 绘制椭球体
% transparency
h.FaceAlpha = transparency;
% 使用固定颜色（例如：红色）
h.FaceColor = RGB; % [0.2, 0.7, 0.2];  % RGB 色彩值，范围是 [0, 1]
% % axis("equal")