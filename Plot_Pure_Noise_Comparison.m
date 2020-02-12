%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference: 'Sparsity-assisted Fault Feature Enhancement: Algorithm-aware versus Model-aware',
% IEEE Transactions on Instrumentation and Measurement, 2020
% Homepage: https://zhaozhibin.github.io/
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2019.6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
clc
clear all
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
rng('default')
rng(19)  

%% Figure initialization
Tstring =  'Time (s)'; 
Fstring = 'Frequency (Hz)';
Astring =  'Amp (m/s^2)';
FontSize = 9;   FontName = 'Times New Roman';
MarkerSize = 4;  LineWidth = 1;
%% Simulation
Fs = 20480;
N = 4096;
Sig_Impulse = QuasiPeiodicImpulseResponse_AM(N, Fs);
t = (0 : N-1) / Fs;
Sig_Cos = 0.0 * cos(2*pi*20*t');

Sigma = 0.6;
Noise = Sigma * randn(N , 1);
Sig_Combine = Sig_Cos + Sig_Impulse' + Noise;




%% Setting the parameters
Q = 2;
r = 5;
J =10;
now = ComputeNow(N,Q,r,J,'radix2');
AH = @(Sig) tqwt_radix2(Sig, Q, r, J);   
A = @(w) itqwt_radix2(w, Q, r , N);
lam = 1.0 * now;
rho = 1;
load Performance_Comparison_Combination_K_Index_Size5_Sigma6.mat
%% method1 : Generalized Structured Shrinkage
[~,Index] = min(mean(GST_Index(:,:,1), 2));
K1 = K(Index); 
Method1.Name = 'WGL';
Method1.Initial_Size = 5;
Method1.SubName = 'MC';
Method1.gamma = 2;
Method1.window = 'gausswin';
z1 = IterGSS(Sig_Combine, A, AH, lam, rho, K1, Method1);


%% method2 : TQWT-L1 based
load Performance_Comparison_Combination_K_Index_Size5_Sigma6.mat
Method2.Name = 'L1';
[~,Index] = min(mean(L1_Index(:,:,1), 2));
K2 = K(Index);   
z2 = IterGSS(Sig_Combine, A, AH, lam, rho, K2, Method2);

%% method3 : Neighbor thresholding
load Performance_Comparison_Combination_K_Index_Size5_Sigma6_NC.mat
params.Q = 2;
params.r = 5;
params.J =10;
[~,Index] = min(mean(NC_Index(:,:,1), 2));
K3 = K(Index);   
z3 = TQWTDe( Sig_Combine, params , 'nc', K3);

%% method4 : Structured Shrinkage
load Performance_Comparison_Combination_K_Index_Size5_Sigma6.mat
Method3.Name = 'WGL';
Method3.Initial_Size = 5;
Method3.SubName = 'L1';
Method3.window = 'gausswin';
[~,Index] = min(mean(ST_Index(:,:,1), 2));
K4 = K(Index);   
z4 = IterGSS(Sig_Combine, A, AH, lam, rho, K4, Method3); 



z1 = real(A(z1));
z2 = real(A(z2));
% z3 = real(A(z3));
z4 = real(A(z4));
%% Calculate RMSE
GST_RMSE = RMSE(z1, Sig_Impulse);
L1_RMSE = RMSE(z2, Sig_Impulse);
NC_RMSE = RMSE(z3', Sig_Impulse);
ST_RMSE = RMSE(z4, Sig_Impulse);


%% Print the time domain
figure(1);clf; 
plot(t, z1, 'b-', 'LineWidth', LineWidth);
title(['GSSA: K=', num2str(K1), ',RMSE=', num2str(round(GST_RMSE*1000)/1000)],'FontSize',FontSize,'FontName',FontName);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Impulse))*1.5;            ylim_max = max(abs(Sig_Impulse))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

xlabel(Tstring);
ylabel(Astring);
filename = ['results', filesep, sprintf('Fig5_a_GSSA.pdf')];
print(filename, '-dpdf');

%% Print the time domain
figure(2);clf;
plot(t, z2, 'b-', 'LineWidth', LineWidth);
title(['BPD: K=', num2str(K2), ',RMSE=', num2str(round(L1_RMSE*1000)/1000)],'FontSize',FontSize,'FontName',FontName);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Impulse))*1.5;            ylim_max = max(abs(Sig_Impulse))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

xlabel(Tstring);
ylabel(Astring);
filename = ['results', filesep, sprintf('Fig5_b_BPD.pdf')];
print(filename, '-dpdf');


%% Print the time domain
figure(3);clf;
plot(t, z4, 'b-', 'LineWidth', LineWidth);
hold off
title(['WGL: K=', num2str(K4), ',RMSE=', num2str(round(ST_RMSE*1000)/1000)],'FontSize',FontSize,'FontName',FontName);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Impulse))*1.5;            ylim_max = max(abs(Sig_Impulse))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

xlabel(Tstring);
ylabel(Astring);

filename = ['results', filesep, sprintf('Fig5_c_WGL.pdf')];
print(filename, '-dpdf');


%% Print the time domain
figure(4);clf;
plot(t, z3, 'b-', 'LineWidth', LineWidth);
title(['NCD: K=', num2str(K3), ',RMSE=', num2str(round(NC_RMSE*1000)/1000)],'FontSize',FontSize,'FontName',FontName);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Impulse))*1.5;            ylim_max = max(abs(Sig_Impulse))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

xlabel(Tstring);
ylabel(Astring);

filename = ['results', filesep, sprintf('Fig5_d_NCD.pdf')];
print(filename, '-dpdf');