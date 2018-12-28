clear all; close all; clc

num_rows = 64;%行数64
num_cols = 64;%列数64
N = 169;      %169个波段
filename = 'hydice';

M = num_rows * num_cols;%每一波段平面，像素个数
fid = fopen([filename '.dat'], 'r','ieee-be');%fid为文件代号，正数表示正确打开，负数表示打开错误，格式为只读，ieee-be
X = fread(fid, [N M],'float32');%读出fid指向的文件M个数据，填入N*M矩阵，浮点数32位
fclose(fid);



%% Global RX算法开始

[N M] = size(X);%将X横纵坐标长度赋给N M
X_mean = mean(X.').'; %按照行求和取平均值（mean本来是按照列求和，转置求均值再转置）
X = X - repmat(X_mean,1,M);%将X_mean扩增[1 M]倍(repmat 复制平铺)，减去均值
Sigma = (X * X')/M; %自相关矩阵 （整体自协方差矩阵）
Sigma_inv = inv(Sigma);  %求逆
for m = 1:M
 D(m) = X(:, m)' * Sigma_inv * X(:, m);  %RX表达式计算式
end
%% Global RX算法结束

r1 = D;%r1为求出的RX算子
tmp = reshape(r1, 64, 64);%构成64*64矩阵
colormap('gray'),imagesc(tmp),colorbar; axis image;%将矩阵数值用颜色表示出来
title('整体法 RX')