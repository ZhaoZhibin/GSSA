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
%% Figure initialization
Tstring1 =  '$\sigma$'; 
Tstring2 =  'K'; 
Astring =  'Average RMSE';
FontSize = 9;   FontName = 'Times New Roman';
MarkerSize = 1;  LineWidth = 1;
%%

load Performance_Comparison_Combination_K_Index_Size5_RMSE.mat
NC = NC_Index;
L1 = L1_Index;
GST = GST_Index;
ST = ST_Index;



%% Print the time domain
figure(1);
hold on
plot(Sigma(1:end), GST(1:end), 'b-*', 'LineWidth', LineWidth);
plot(Sigma(1:end), L1(1:end) , 'r-^', 'LineWidth', LineWidth);
plot(Sigma(1:end), ST(1:end) , 'g-o', 'LineWidth', LineWidth);
plot(Sigma(1:end), NC(1:end) , 'k->', 'LineWidth', LineWidth);
hold off
box on
legend1 = legend('GSSA' , 'BPD', 'WGL', 'NCD');
set(legend1,'location','best','Orientation','horizontal', 'FontSize',FontSize,'FontName',FontName)
legend boxoff
xlim_min = 0.2; xlim_max = 0.6;
ylim_min = min(abs(GST))*0.8;            ylim_max = max(abs(NC)) * 1.2;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

xlabel(Tstring1);
ylabel(Astring);
filename = ['results', filesep, sprintf('Fig4_b.pdf')];
print(filename, '-dpdf');

%%
load Performance_Comparison_Combination_K_Index_Size5_Sigma6_NC.mat
load Performance_Comparison_Combination_K_Index_Size5_Sigma6.mat
NC = mean(NC_Index(1:100,:,1), 2);
L1 = mean(L1_Index(:,:,1), 2);
GST = mean(GST_Index(:,:,1), 2);
ST = mean(ST_Index(:,:,1), 2);

% Print the time domain
figure(2);
hold on
plot(K(1:4:end), GST(1:4:end), 'b-*', 'LineWidth', LineWidth);
plot(K(1:4:end), L1(1:4:end),'r-^', 'LineWidth', LineWidth);
plot(K(1:4:end), ST(1:4:end),'g-o', 'LineWidth', LineWidth);
plot(K(1:4:end), NC(1:4:end),'k->', 'LineWidth', LineWidth);
hold off
box on
legend1 = legend('GSSA' , 'BPD', 'WGL', 'NCD');
set(legend1,'location','best','Orientation','horizontal', 'FontSize',FontSize,'FontName',FontName)
legend boxoff
xlim_min = 0; xlim_max = max(K);
ylim_min = min(abs(L1))*0.8;            ylim_max = max(abs(GST)) * 1.1;  
xylim = [xlim_min,xlim_max,ylim_min,ylim_max]; axis(xylim);

xlabel(Tstring2);
ylabel(Astring);
filename = ['results', filesep, sprintf('Fig4_a.pdf')];
print(filename, '-dpdf');
