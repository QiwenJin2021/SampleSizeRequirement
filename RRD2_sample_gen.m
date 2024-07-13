
clear all;clc;
addpath(genpath(pwd));
filepath = 'D:\Mycode\ParticleAnalysis\PSD_simu_v2\SampleData\RRD2\';

type = 'RRD2';
method = 'ARM';
n = 1.2:0.2:6;
xmin = 1;
xmax = 101;
De = 50;
P = 0.95;

% Generation of GSD2 samples
N = 10000000;

for j = 1:length(n)
    fprintf('n=%2d\n', n(j));
    X = gen_sample(method,type,N,n(j),xmin,xmax,De);
    fname = [filepath 'RRD2_n_'  num2str(n(j)) '.mat'];
    save(fname, 'X' );

end

