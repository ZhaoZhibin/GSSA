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
rng(45)
%% Figure initialization
FontSize = 10;   FontName = 'Times New Roman';
MarkerSize = 4;  LineWidth = 1;
%%
FlagFigureAutoSave = 1;
currentFolder = pwd;

%% Penalty
%% Generate the data
x = -5 : 0.01 : 5;
lambda = 1;

%% six different penalty
% L1
phi1 = penalty( 'L1', lambda );
y1 = phi1(x);
% L0
phi2 = penalty( 'L0', lambda );
y2 = phi2(x);
% Lp(p=0.5)
phi3 = penalty( 'Lp', lambda );
y3 = phi3(x, 0.5);
% SCAD
phi4 = penalty( 'SCAD', lambda );
y4 = phi4(x, 3.7);
% MC
phi5 = penalty( 'MC', lambda );
y5 = phi5(x, 2);
%% Print the time domain
figure();clf; 
hold on
ph(1) = plot(x, y1, 'b-d', 'LineWidth', LineWidth,'MarkerIndices',750);
ph(2) = plot(x, y2, 'r-s', 'LineWidth', LineWidth,'MarkerIndices',630);
ph(3) = plot(x, y3, 'g->', 'LineWidth', LineWidth,'MarkerIndices',750);
ph(4) = plot(x, y4, 'k-*', 'LineWidth', LineWidth,'MarkerIndices',750);
ph(5) = plot(x, y5, 'm-o', 'LineWidth', LineWidth,'MarkerIndices',600);
hold off
box on
legend1 = legend(ph, '$l_1$' , '$l_0$' , '$l_{1/2}$', 'SCAD (a=3.7)', 'MC ($\gamma=2$)');
set(legend1,'location','north','Orientation','vertical', 'FontSize',FontSize+2,'FontName',FontName,'Interpreter','latex')
legend boxoff

xlim_min = min(-5); xlim_max = max(x);
ylim_min = 0;            ylim_max = 5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

filename = ['results', filesep, sprintf('Demo1_Penalties.pdf')];
print(filename, '-dpdf');

%% Shrinkage
%% six different Shrinkage
% L1
params.lambda = lambda;
threshold1 = Shrinkage( 'L1', params );
s1 = threshold1(x);
% L0
params.lambda = lambda^2/2;
threshold2 = Shrinkage( 'L0', params );
s2 = threshold2(x);
% Lp(p=1/2)
params.lambda = sqrt((4 / (54)^(1/3))^3);
params.p = 1/2;
threshold3 = Shrinkage( 'Lp', params );
s3 = threshold3(x);
% SCAD
params.lambda = lambda;
params.a = 3.7;
threshold4 = Shrinkage( 'SCAD', params );
s4 = threshold4(x);
% MC
params.lambda = lambda;
params.gamma = 2;
threshold5 = Shrinkage( 'MC', params );
s5 = threshold5(x);
%% Print the time domain
figure();
hold on
ph(1) = plot(x, s1, 'b-d', 'LineWidth', LineWidth,'MarkerIndices',750);
ph(2) = plot(x, s2, 'r-s', 'LineWidth', LineWidth,'MarkerIndices',630);
ph(3) = plot(x, s3, 'g->', 'LineWidth', LineWidth,'MarkerIndices',750);
ph(4) = plot(x, s4, 'k-*', 'LineWidth', LineWidth,'MarkerIndices',750);
ph(5) = plot(x, s5, 'm-o', 'LineWidth', LineWidth,'MarkerIndices',650);
hold off
box on
legend1 = legend(ph, '$l_1$' , '$l_0$' , '$l_{1/2}$', 'SCAD (a=3.7)', 'MC ($\gamma=2$)');
set(legend1,'location','best','Orientation','vertical', 'FontSize',FontSize+2,'FontName',FontName,'Interpreter','latex')
legend boxoff

xlim_min = min(0); xlim_max = max(x);
ylim_min = 0;            ylim_max = 5;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

filename = ['results', filesep, sprintf('Demo1_Shrinkage.pdf')];
print(filename, '-dpdf');