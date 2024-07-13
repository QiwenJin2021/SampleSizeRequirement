
% Gy采样理论计算
% std = C*d^3*(1/MS-1/ML)
%

clear all;clc;
addpath(genpath(pwd));

type = 'RRD2';
method = 'ARM';
n = 1.2:0.2:6;
xmin = 1;
xmax = 101;
De = 50;
P = 0.95;
% filepath = 'D:\Mycode\ParticleAnalysis\PSD_simu_v2\SampleData\GSD2_ARM\';
filepath = ['D:\Mycode\ParticleAnalysis\PSD_simu_v2\SampleData\', type, '_', method,  '\'];

% Factors in Gy theory
% FSE = fgbc*d^3(1/Ms-1/ML)
% f:shape factor, 0.5 for spherical particle
% g:size distribution factor, 0.25 for wide-size distribution
% b:liberation factor, 1 for particles are completely liberated
% c:constitution factor, c=(1-aL/a)^2/(aL/a)*rhoc+(1-aL/a)rhom
% aL:average concentration of the lot, 10%
% a:concentration of the analyte in the critical particles, 100%
% rhoc:density of the critical particles, 1.5g/cm^3
% rhom: density of the matrix, 1.5g/cm^3
f = 0.5;
g = 0.25;
b = 1;
aL = 0.1;
a = 1;
rhoc = 1.5; %g/cm^3
rhom = 1.5; %g/cm^3
c = (1-aL/a)^2/(aL/a)*rhoc + (1-aL/a)*rhom; % g/cm^3
C = f*g*b*c; % g/cm^3

times = 2000;


% 固定样本量，改变n
SampleSize = 100000;
SampleSize = SampleSize';
len_SampleSize = length(SampleSize);
for kk = 1:length(n)
    %     [xm,xn] = get_bin(type,n(kk),xmin,xmax,De);
    %     D_the = calc_D_theory(type,n(kk),xmin,xmax,De);
    Dm95 = 0.95^(1/n(kk))*(xmax-xmin)+xmin; %um
    fname = [filepath type '_n_'  num2str(n(kk)) '.mat'];
    load(fname);
    ML = sum(X.^3)*4/3*pi*rhoc*10e-12; %g

    for j = 1:times
        X1 = datasample(X,SampleSize);
        MS = sum(X1.^3)*4/3*pi*rhoc*10e-12;
        mass_per(j,kk) = MS/ML;
        std_gy(j,kk) = C*Dm95^3*10e-12*(1/MS-1/ML);
    end
    std_gy_avg(kk,1) = mean(std_gy(:,kk));
    std_gy_std(kk,1) = std(std_gy(:,kk));
end

% 固定n，改变样本量
% n = 6;
% Dm95 = 0.95^(1/n)*(xmax-xmin)+xmin;
% fname = [filepath type '_n_'  num2str(n) '.mat'];
% load(fname);
% ML = sum(X.^3)*4/3*pi*rhoc*10e-12;
% 
% SampleSize = [1000:1000:20000,30000:10000:100000];
% SampleSize = SampleSize';
% 
% for i = 1:length(SampleSize)
%     fprintf('i=%2d\n',i);
%     for j = 1:times
%         X1 = datasample(X,SampleSize(i));
%         MS = sum(X1.^3)*4/3*pi*rhoc*10e-12;
%         mass_per(j,i) = MS/ML;
%         std_gy(j,i) = C*Dm95^3*10e-12*(1/MS-1/ML);
%     end
%     std_gy_avg(i,1) = mean(std_gy(:,i));
%     std_gy_std(i,1) = std(std_gy(:,i));
% end



