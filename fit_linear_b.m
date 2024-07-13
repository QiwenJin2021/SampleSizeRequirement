clear all;clc;
addpath(genpath(pwd));

type = 'GSD2';
method = 'ARM';
n = 1.2:0.2:6;
xmin = 1;
xmax = 101;
De = 50;
P = 0.95;

for kk = 1:length(n)
    
    SR(kk,:) = calc_SR(type,n(kk),xmin,xmax,De);

end
SRm = log(SR(:,1));

filepath = ['D:\Mycode\ParticleAnalysis\PSD_simu_v2\SampleData\' type '_' method '\'];

matfile = [filepath 'Results2_' type '.mat'];
load(matfile);

for i = 2:5
    b = B(:,i);
    [p,S] = polyfit(SRm,b,6);
%     powerFitType = fittype('k0+k4*SRm^4+k3*SRm^3+k2*SRm^2+k1*SRm',...
%         'dependent',{'b'},'independent',{'SRm'},...
%         'coefficients',{'k0','k4','k3','k2','k1',});
%     powerFitType = fittype('k0+k2*SRm^2+k1*SRm',...
%         'dependent',{'b'},'independent',{'SRm'},...
%         'coefficients',{'k0','k2','k1',});
%     [power_fit,gof] = fit(b,SRm,powerFitType);
    %         A(j,i) = power_fit.a;
%     B(j,i) = power_fit.b;
%     rs2(j,i) = gof.rsquare;
end