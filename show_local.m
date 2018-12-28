clear all; close all; clc

num_rows = 64;%����64
num_cols = 64;%����64
N = 169;      
filename = 'hydice';

M = num_rows * num_cols;%ÿһ����ƽ�棬���ظ���
fid = fopen([filename '.dat'], 'r','ieee-be');%fidΪ�ļ����ţ�������ʾ��ȷ�򿪣�������ʾ�򿪴���
X = fread(fid, [N M],'float32');%����fidָ����ļ�M�����ݣ�����N*M���󣬸�����32λ
fclose(fid);

X_mean = mean(X.').'; %���������ȡƽ��ֵ��mean�����ǰ�������ͣ�ת�����ֵ��ת�ã�
X = X - repmat(X_mean,1,M);%��X_mean����[1 M]��(repmat ����ƽ��)����ȥ��ֵ
DataTest = reshape(X', num_rows, num_cols, N);%ת�ã����г���ά64*64*169
DataTest_ori = DataTest(:, :, :);

win_out = 11; % �ⴰ��
win_in = 3; % �ڴ���


%% Local RX�㷨��ʼ

Data = DataTest_ori;
[a b c] = size(Data);%�������ݼ��м���
result = zeros(a, b);%����һ��a��,b�е������
t = fix(win_out/2);%���㷽��ȡ��,Z_OUT
t1 = fix(win_in/2);%Z_IN
M1 = win_out^2;
% padding avoid edges�������˿��Ե��������ħ����һ�棩�����Ҫ��Ϊ����
DataTest = zeros(3*a, 3*b, c);
DataTest(a+1:2*a, b+1:2*b, :) = Data;%ԭʼ��������м�
DataTest(a+1:2*a, 1:b, :) = Data(:, b:-1:1, :);%ԭʼ�������ҷ�ת��������
DataTest(a+1:2*a, 2*b+1:3*b, :) = Data(:, b:-1:1, :);%ԭʼ�������ҷ�ת�����붫
DataTest(1:a, :, :) = DataTest(2*a:-1:(a+1), :, :);%�м���о��������������
DataTest(2*a+1:3*a, :, :) = DataTest(2*a:-1:(a+1), :, :);%�м���о��������������

for i = 1+b: 2*b 
    for j = 1+a: 2*a
        block = DataTest(j-t: j+t, i-t: i+t, :);%�涨�ⴰ�ڴ�С(11*11)*169
        y = squeeze(DataTest(j, i, :)).';%��1*1*169����ά����ת��1*169�Ķ�ά����y���ε���ÿ�����ص�Ķ�ά����squeeze������
        block(t-t1+1:t+t1+1, t-t1+1:t+t1+1, :) = NaN;%���ڴ��ڴ���3*3��*169�����ڸ�Ϊ����ֵԪ��NaN
        block = reshape(block, M1, c);%����άתΪ��ά121*169��ÿһ��Ϊһ������
        block(isnan(block(:, 1)), :) = [];%���վ��󸳸�nanֵ��ȥ��ÿһ�е�nanֵ������ÿһ��Ϊ121-9=112������(isnan�ж��Ƿ�Ϊ����ֵԪ�أ�����Ƿ���1)
        H = block.';  % ת�ã�ÿһ��Ϊһ������
        Sigma = (H * H');%ͬglobal�ļ��㣬����ؾ���
        Sigma_pinv = pinv(Sigma); %ͬglobal�ļ��㣬���������
%         Sigma_inv = inv(Sigma);
       result(j-a, i-b) = y * Sigma_pinv * y';%Rx����
%         result(j-a, i-b) = y * Sigma_inv * y';
    end
end
%% Local RX�㷨�������ó�RX����result
r2 = result;
r2 = reshape(r2, 1, M);%��ÿ�����ص��Rx���ֵ�ų��о���
tmp = reshape(r2, 64, 64);%�����64*64
colormap('gray'),imagesc(tmp),colorbar; axis image;%����ֵ�Ĵ�С��ɫ��ͬ��ͼ
title('˫���� RX')
