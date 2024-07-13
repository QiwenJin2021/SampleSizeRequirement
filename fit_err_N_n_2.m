
clear all;clc;
addpath(genpath(pwd));

type = 'RRD2';
method = 'ARM';
n = 1.2:0.2:6;
xmin = 1;
xmax = 101;
De = 50;
P = 0.95;

filepath = ['D:\Mycode\ParticleAnalysis\PSD_simu_v2\SampleData\' type '_' method '\'];

N = [1000:1000:20000,30000:10000:100000];
N_log = log10(N);
N_log = N_log';

% for i = 1:length(N)
%     for j = 1:length(n)
%         fprintf('n=%2d\n', n(j));
%         matfile = [filepath 'Results_' type '_n_' num2str(n(j)) '.mat'];
%         load(matfile);
%         err_P_N(j,:) = err_P(i,:);
%         err_avg_N(j,:) = err_avg(i,:);
% 
%     end
%     file = [filepath 'Results_' type '_N_' num2str(N(i)) '.mat'];
%     save(file,'err_P_N','err_avg_N');
% end

for j = 1:length(n)
    fprintf('n=%2d\n', n(j));
    matfile = [filepath 'Results_' type '_n_' num2str(n(j)) '.mat'];
    load(matfile);
    for i = 1:5
        err_P_log = log10(err_P(:,i));
        err_avg_log = log10(err_avg(:,i));
        powerFitType = fittype('-2 * err_P_log + b',...
            'dependent',{'N_log'},'independent',{'err_P_log'},...
            'coefficients',{'b'});
        [power_fit,gof] = fit(err_P_log,N_log,powerFitType);
%         A(j,i) = power_fit.a;
        B(j,i) = power_fit.b;
        rs2(j,i) = gof.rsquare;
        rs2_adj(j,i) = gof.adjrsquare;
    end

end
file = [filepath 'FitResults_k=-2_' type '.mat'];
save(file,'B','rs2','rs2_adj');

