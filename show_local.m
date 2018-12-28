clear all; close all; clc

num_rows = 64;%行数64
num_cols = 64;%列数64
N = 169;      
filename = 'hydice';

M = num_rows * num_cols;%每一波段平面，像素个数
fid = fopen([filename '.dat'], 'r','ieee-be');%fid为文件代号，正数表示正确打开，负数表示打开错误
X = fread(fid, [N M],'float32');%读出fid指向的文件M个数据，填入N*M矩阵，浮点数32位
fclose(fid);

X_mean = mean(X.').'; %按照行求和取平均值（mean本来是按照列求和，转置求均值再转置）
X = X - repmat(X_mean,1,M);%将X_mean扩增[1 M]倍(repmat 复制平铺)，减去均值
DataTest = reshape(X', num_rows, num_cols, N);%转置，排列成三维64*64*169
DataTest_ori = DataTest(:, :, :);

win_out = 11; % 外窗口
win_in = 3; % 内窗口


%% Local RX算法开始

Data = DataTest_ori;
[a b c] = size(Data);%返回数据几行几列
result = zeros(a, b);%创建一个a行,b列的零矩阵
t = fix(win_out/2);%朝零方向取整,Z_OUT
t1 = fix(win_in/2);%Z_IN
M1 = win_out^2;
% padding avoid edges，（填充八块边缘，想三阶魔方的一面），填充要求互为镜像
DataTest = zeros(3*a, 3*b, c);
DataTest(a+1:2*a, b+1:2*b, :) = Data;%原始数据填充中间
DataTest(a+1:2*a, 1:b, :) = Data(:, b:-1:1, :);%原始数据左右翻转，填入西
DataTest(a+1:2*a, 2*b+1:3*b, :) = Data(:, b:-1:1, :);%原始数据左右翻转，填入东
DataTest(1:a, :, :) = DataTest(2*a:-1:(a+1), :, :);%中间的行镜像填入上面的行
DataTest(2*a+1:3*a, :, :) = DataTest(2*a:-1:(a+1), :, :);%中间的行镜像填入下面的行

for i = 1+b: 2*b 
    for j = 1+a: 2*a
        block = DataTest(j-t: j+t, i-t: i+t, :);%规定外窗口大小(11*11)*169
        y = squeeze(DataTest(j, i, :)).';%将1*1*169的三维矩阵转成1*169的二维矩阵。y依次等于每个像素点的二维矩阵（squeeze函数）
        block(t-t1+1:t+t1+1, t-t1+1:t+t1+1, :) = NaN;%将内窗口处（3*3）*169区域内赋为非数值元素NaN
        block = reshape(block, M1, c);%将三维转为二维121*169，每一列为一个波段
        block(isnan(block(:, 1)), :) = [];%将空矩阵赋给nan值（去掉每一列的nan值，现在每一列为121-9=112个数）(isnan判断是否为非数值元素，如果是返回1)
        H = block.';  % 转置，每一行为一个波段
        Sigma = (H * H');%同global的计算，自相关矩阵
        Sigma_pinv = pinv(Sigma); %同global的计算，计算逆矩阵
%         Sigma_inv = inv(Sigma);
       result(j-a, i-b) = y * Sigma_pinv * y';%Rx算子
%         result(j-a, i-b) = y * Sigma_inv * y';
    end
end
%% Local RX算法结束，得出RX算子result
r2 = result;
r2 = reshape(r2, 1, M);%将每个像素点的Rx结果值排成行矩阵
tmp = reshape(r2, 64, 64);%重组成64*64
colormap('gray'),imagesc(tmp),colorbar; axis image;%根据值的大小颜色不同画图
title('双窗口 RX')
