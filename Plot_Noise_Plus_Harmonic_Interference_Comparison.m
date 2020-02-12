%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference: 'Sparsity-assisted Fault Feature Enhancement: Algorithm-aware versus Model-aware',
% IEEE Transactions on Instrumentation and Measurement, 2020
% Homepage: https://zhaozhibin.github.io/
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2019.6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
A2string =  'Amp (m^2/s^4)';
FontSize = 9;   FontName = 'Times New Roman';
MarkerSize = 3;  LineWidth = 1;
%% Simulation
Fs = 20480;
N = 4096;
Sig_Impulse = QuasiPeiodicImpulseResponse_AM(N, Fs);
t = (0 : N-1) / Fs;
Sig_Cos = 0.5 * cos(2*pi*160*t') + 0.3 * cos(2*pi*320*t');

Sigma = 0.6;
Noise = Sigma * randn(N , 1);
Sig_Combine = Sig_Cos + Sig_Impulse' + Noise;
% Sig_Combine = Sig_Cos + Sig_Impulse + Noise;
[ yf2, f2 ] = Dofft( Sig_Combine , Fs , 0);
[ yfh, fh ] = Hilbert_envelope( Sig_Combine , Fs , 1);


%% Setting the parameters
Q = 2;
r = 5;
J =10;
now = ComputeNow(N,Q,r,J,'radix2');
AH = @(Sig) tqwt_radix2(Sig, Q, r, J);   
A = @(w) itqwt_radix2(w, Q, r , N);
lam = 1.0 * now;
rho = 1;
load Harmonic_Inference_Best_K.mat
%% method1 : Generalized Structured Shrinkage
[~,Index] = min(GST_RMSE);
K1 = K(Index); 
Method1.Name = 'WGL';
Method1.Initial_Size = 5;
Method1.SubName = 'MC';
Method1.gamma = 2;
Method1.window = 'gausswin';
z1 = IterGSS(Sig_Combine, A, AH, lam, rho, K1, Method1);


%% method2 : TQWT-L1 based
Method2.Name = 'L1';
[~,Index] = min(L1_RMSE);
K2 = K(Index);   
z2 = IterGSS(Sig_Combine, A, AH, lam, rho, K2, Method2);

%% method3 : Neighbor thresholding
params.Q = 2;
params.r = 5;
params.J =10;
[~,Index] = min(NC_RMSE);
K3 = K(Index);   
z3 = TQWTDe( Sig_Combine, params , 'nc', K3);

%% method4 : Structured Shrinkage
Method3.Name = 'WGL';
Method3.Initial_Size = 5;
Method3.SubName = 'L1';
Method3.window = 'gausswin';
[~,Index] = min(ST_RMSE);
K4 = K(Index);   
z4 = IterGSS(Sig_Combine, A, AH, lam, rho, K4, Method3); 


%% Time domain
z1 = real(A(z1));
z2 = real(A(z2));
z4 = real(A(z4));

%% Frequency domain
[ yfz1, fz1 ] = Hilbert_envelope( z1 , Fs , 1);
[ yfz2, fz2 ] = Hilbert_envelope( z2 , Fs , 1);
[ yfz3, fz3 ] = Hilbert_envelope( z3 , Fs , 1);
[ yfz4, fz4 ] = Hilbert_envelope( z4 , Fs , 1);


%% Calculate RMSE
GST_RMSE = RMSE(z1, Sig_Impulse);
L1_RMSE = RMSE(z2, Sig_Impulse);
NC_RMSE = RMSE(z3', Sig_Impulse);
ST_RMSE = RMSE(z4, Sig_Impulse);

%% Print the time domain
figure;
subplot(221)
plot(t, Sig_Impulse, 'b-', 'LineWidth', LineWidth);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Impulse))*1.5;            ylim_max = max(abs(Sig_Impulse))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Tstring);
ylabel(Astring);
subplot(222)
plot(t, Sig_Cos, 'b-', 'LineWidth', LineWidth);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Cos))*1.5;            ylim_max = max(abs(Sig_Cos))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Tstring);
ylabel(Astring);
subplot(223)
plot(t, Sig_Combine, 'b-', 'LineWidth', LineWidth);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Combine))*1.5;            ylim_max = max(abs(Sig_Combine))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Tstring);
ylabel(Astring);
subplot(224)
plot(f2, yf2, 'b-', 'LineWidth', LineWidth);
xlim_min = 0; xlim_max = 5000;
ylim_min = 0;            ylim_max = max(abs(yf2))*1.3;   
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Fstring);
ylabel(Astring);
filename = ['results', filesep, sprintf('Fig6_Simulation.pdf')];
print(filename, '-dpdf');




%% Print GSSA
figure;
subplot(121)
plot(t, z1, 'b-', 'LineWidth', LineWidth);
title(['GSSA: K=', num2str(K1), ',RMSE=', num2str(round(GST_RMSE*1000)/1000)],'FontSize',FontSize,'FontName',FontName);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Impulse))*1.5;            ylim_max = max(abs(Sig_Impulse))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Tstring);
ylabel(Astring);
subplot(122)
plot(fz1, yfz1, 'b-', 'LineWidth', LineWidth);
hold on
% impulsive frequencies
plot(100, yfz1(fz1==100), 'ro', 'MarkerSize', MarkerSize);
plot(200, yfz1(fz1==200), 'ro', 'MarkerSize', MarkerSize);
plot(300, yfz1(fz1==300), 'ro', 'MarkerSize', MarkerSize);
% harmonic frequencies
plot(160, yfz1(fz1==160), 'g^', 'MarkerSize', MarkerSize);
plot(320, yfz1(fz1==320), 'g^', 'MarkerSize', MarkerSize);
plot(480, yfz1(fz1==480), 'g^', 'MarkerSize', MarkerSize);
hold off
xlim_min = 0; xlim_max = 1000;
ylim_min = 0;            ylim_max = max(abs(yfz1))*1.3;   
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Fstring);
ylabel(A2string);
filename = ['results', filesep, sprintf('Fig7_ab_GSSA.pdf')];
print(filename, '-dpdf');


%% Print BPD

figure;
subplot(121)
plot(t, z2, 'b-', 'LineWidth', LineWidth);
title(['BPD: K=', num2str(K2), ',RMSE=', num2str(round(L1_RMSE*1000)/1000)],'FontSize',FontSize,'FontName',FontName);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Impulse))*1.5;            ylim_max = max(abs(Sig_Impulse))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Tstring);
ylabel(Astring);
subplot(122)
plot(fz2, yfz2, 'b-', 'LineWidth', LineWidth);
hold on
% impulsive frequencies
plot(100, yfz2(fz2==100), 'ro', 'MarkerSize', MarkerSize);
plot(200, yfz2(fz2==200), 'ro', 'MarkerSize', MarkerSize);
plot(300, yfz2(fz2==300), 'ro', 'MarkerSize', MarkerSize);
% harmonic frequencies
plot(160, yfz2(fz2==160), 'g^', 'MarkerSize', MarkerSize);
plot(320, yfz2(fz2==320), 'g^', 'MarkerSize', MarkerSize);
plot(480, yfz2(fz2==480), 'g^', 'MarkerSize', MarkerSize);
hold off
xlim_min = 0; xlim_max = 1000;
ylim_min = 0;            ylim_max = max(abs(yfz2))*1.3;   
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Fstring);
ylabel(A2string);
filename = ['results', filesep, sprintf('Fig7_cd_BPD.pdf')];
print(filename, '-dpdf');


%% Print WGL

figure;
subplot(121)
plot(t, z4, 'b-', 'LineWidth', LineWidth);
title(['WGL: K=', num2str(K4), ',RMSE=', num2str(round(ST_RMSE*1000)/1000),'0'],'FontSize',FontSize,'FontName',FontName);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Impulse))*1.5;            ylim_max = max(abs(Sig_Impulse))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Tstring);
ylabel(Astring);
subplot(122)
%% Print the square envelope spretrum
plot(fz4, yfz4, 'b-', 'LineWidth', LineWidth);
hold on
% impulsive frequencies
plot(100, yfz4(fz4==100), 'ro', 'MarkerSize', MarkerSize);
plot(200, yfz4(fz4==200), 'ro', 'MarkerSize', MarkerSize);
plot(300, yfz4(fz4==300), 'ro', 'MarkerSize', MarkerSize);
% harmonic frequencies
plot(160, yfz4(fz4==160), 'g^', 'MarkerSize', MarkerSize);
plot(320, yfz4(fz4==320), 'g^', 'MarkerSize', MarkerSize);
plot(480, yfz4(fz4==480), 'g^', 'MarkerSize', MarkerSize);
hold off
xlim_min = 0; xlim_max = 1000;
ylim_min = 0;            ylim_max = max(abs(yfz4))*1.3;   
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Fstring);
ylabel(A2string);
filename = ['results', filesep, sprintf('Fig7_ef_WGL.pdf')];
print(filename, '-dpdf');


%% Print NCD
figure;
subplot(121)
plot(t, z3, 'b-', 'LineWidth', LineWidth);
title(['NCD: K=', num2str(K3), ',RMSE=', num2str(round(NC_RMSE*1000)/1000)],'FontSize',FontSize,'FontName',FontName);
xlim_min = min(t); xlim_max = max(t);
ylim_min = -max(abs(Sig_Impulse))*1.5;            ylim_max = max(abs(Sig_Impulse))*1.5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Tstring);
ylabel(Astring);
subplot(122)
%% Print the square envelope spretrum
plot(fz3, yfz3, 'b-', 'LineWidth', LineWidth);
hold on
% impulsive frequencies
plot(100, yfz3(fz3==100), 'ro', 'MarkerSize', MarkerSize);
plot(200, yfz3(fz3==200), 'ro', 'MarkerSize', MarkerSize);
plot(300, yfz3(fz3==300), 'ro', 'MarkerSize', MarkerSize);
% harmonic frequencies
plot(160, yfz3(fz3==160), 'g^', 'MarkerSize', MarkerSize);
plot(320, yfz3(fz3==320), 'g^', 'MarkerSize', MarkerSize);
plot(480, yfz3(fz3==480), 'g^', 'MarkerSize', MarkerSize);
hold off
xlim_min = 0; xlim_max = 1000;
ylim_min = 0;            ylim_max = max(abs(yfz3))*1.3;   
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);
xlabel(Fstring);
ylabel(A2string);
filename = ['results', filesep, sprintf('Fig7_gh_NCD.pdf')];
print(filename, '-dpdf');
