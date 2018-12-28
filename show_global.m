clear all; close all; clc

num_rows = 64;%����64
num_cols = 64;%����64
N = 169;      %169������
filename = 'hydice';

M = num_rows * num_cols;%ÿһ����ƽ�棬���ظ���
fid = fopen([filename '.dat'], 'r','ieee-be');%fidΪ�ļ����ţ�������ʾ��ȷ�򿪣�������ʾ�򿪴��󣬸�ʽΪֻ����ieee-be
X = fread(fid, [N M],'float32');%����fidָ����ļ�M�����ݣ�����N*M���󣬸�����32λ
fclose(fid);



%% Global RX�㷨��ʼ

[N M] = size(X);%��X�������곤�ȸ���N M
X_mean = mean(X.').'; %���������ȡƽ��ֵ��mean�����ǰ�������ͣ�ת�����ֵ��ת�ã�
X = X - repmat(X_mean,1,M);%��X_mean����[1 M]��(repmat ����ƽ��)����ȥ��ֵ
Sigma = (X * X')/M; %����ؾ��� ��������Э�������
Sigma_inv = inv(Sigma);  %����
for m = 1:M
 D(m) = X(:, m)' * Sigma_inv * X(:, m);  %RX���ʽ����ʽ
end
%% Global RX�㷨����

r1 = D;%r1Ϊ�����RX����
tmp = reshape(r1, 64, 64);%����64*64����
colormap('gray'),imagesc(tmp),colorbar; axis image;%��������ֵ����ɫ��ʾ����
title('���巨 RX')