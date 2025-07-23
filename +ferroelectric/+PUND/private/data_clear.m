function [x_clean,y_clean,z_clean] = data_clear(x,y,z)
%   剔除掉变化率特别小的x
%   此处显示详细说明
dx = diff(x);
dy = diff(y);

% 计算斜率（处理 dx 非零的情况）
slope = dy ./ dx;

% 设置阈值（避免 dx 太小造成除以零）
dx_thresh = 0.01;     % x 变化非常小
% 找出异常点：x 几乎不变，但 y 变化大（即斜率极大）
is_outlier = abs(dx) < dx_thresh;

% 第四步：剔除异常点
x_clean = x(~is_outlier);
y_clean = y(~is_outlier);
z_clean = z(~is_outlier);
end
