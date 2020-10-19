%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
clc
clear all
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
rng('default')
rng(19)  

%% Figure initialization
global FontSize FontName;
Tstring =  'Time (s)'; 
Fstring = 'Frequency (Hz)';
Astring =  'Amp (m/s^2)';
FontSize = 11;   FontName = 'Times New Roman';
MarkerSize = 4;  LineWidth = 1;
%%
FlagFigureAutoSave = 1;
currentFolder = pwd;
%% Simulation
Fs = 20480;
N = 4096;
mode = 'outer';
% [Sig_Impulse , t] = MakeSignalBearing( Fs, N, 'outer');
Sig_Impulse = QuasiPeiodicImpulseResponse_AM(N, Fs);
t = (0 : N-1) / Fs;
Sig_Cos = 0.5 * cos(2*pi*160*t') + 0.3 * cos(2*pi*320*t');

Sigma = 0.6;
Noise = Sigma * randn(N , 1);
Sig_Combine = Sig_Cos + Sig_Impulse' + Noise;
% Sig_Combine = Sig_Cos + Sig_Impulse + Noise;
[ yf2, f2 ] = Dofft( Sig_Combine , Fs , 0);

K = 10 : 10 : 1000;
for i = 1 : length(K)
    %% Setting the parameters
    i
    Q = 2;
    r = 5;
    J =10;
    now = ComputeNow(N,Q,r,J,'radix2');
    AH = @(Sig) tqwt_radix2(Sig, Q, r, J);   
    A = @(w) itqwt_radix2(w, Q, r , N);
    lam = 1.0 * now;
    rho = 1;
    %% method1 : Generalized Structured Shrinkage
    K1 = K(i); 
    Method1.Name = 'WGL';
    Method1.Initial_Size = 5;
    Method1.SubName = 'MC';
    Method1.gamma = 2;
    Method1.window = 'gausswin';
    z1 = IterGSS(Sig_Combine, A, AH, lam, rho, K1, Method1);


    %% method2 : TQWT-L1 based
    Method2.Name = 'L1';
    K2 = K(i);   
    z2 = IterGSS(Sig_Combine, A, AH, lam, rho, K2, Method2);

    %% method3 : Neighbor thresholding
    params.Q = 2;
    params.r = 5;
    params.J =10;
    K3 = K(i);   
    z3 = TQWTDe( Sig_Combine, params , 'nc', K3);

    %% method4 : Structured Shrinkage
    Method3.Name = 'WGL';
    Method3.Initial_Size = 5;
    Method3.SubName = 'L1';
    Method3.window = 'gausswin';
    K4 = K(i);   
    z4 = IterGSS(Sig_Combine, A, AH, lam, rho, K4, Method3); 
    z1 = real(A(z1));
    z2 = real(A(z2));
    % z3 = real(A(z3));
    z4 = real(A(z4));
    %% Calculate RMSE
    GST_RMSE(i) = RMSE(z1, Sig_Impulse);
    L1_RMSE(i) = RMSE(z2, Sig_Impulse);
    NC_RMSE(i) = RMSE(z3', Sig_Impulse);
    ST_RMSE(i) = RMSE(z4, Sig_Impulse);
end
save('Harmonic_Inference_Best_K.mat', 'K', 'GST_RMSE', 'L1_RMSE', 'NC_RMSE', 'ST_RMSE')




