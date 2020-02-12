function [ DSig ] = TQWTDe( Sig, params, Method, K)
% The function realizes that TQWT threshold denoising
% Input:
%         Sig : the original signal
%         params : the parameters of TQWT
%         Method : the threshold method('h':hardthreshold; 's':soft-threshold; 'nc':neighbor coefficients threshold)
%         Sigma : the standard deviation of the signal
% Output:
%         DSig : the signal after wavelet denoising
addpath tqwt_matlab_toolbox
if ~exist('K', 'var') || isempty(K)
    K = 0.02;                           
end
Q = params.Q;
r = params.r;
J = params.J;
Sig = Sig(:);
w = tqwt_radix2(Sig, Q, r, J);                                          % Wavelet decomposition
Temp = [];
normA = ComputeNow(length(Sig),Q,r,J,'radix2');
for i = 1:numel(w)
    Temp = [Temp ; w{i}(:) / normA(i)];
end    
T_Value = K_sparsity(Temp, K); 
for i = 1: J+1                                                       % Denoising scale by scale
    wi = w{i};    
    % Calculate the threshold according to the work of Donoho
    tau = T_Value * normA(i);
    if strcmp(Method , 's') || strcmp(Method, 'h')
    % Using the function inside Matlab to denosing(We provide the codes below)
    coefi = wthresh(wi, Method , tau);
    else
    coefi = wth_nc(wi , T_Value*normA(i))';
    end
    w{i} = coefi;
end

DSig = itqwt_radix2(w, Q, r , length(Sig))';
end

function [Dcoef] = wth_nc(coef , T)
% The function realizes the neighbor coefficients threshold denoising
% Input:
%          coef : the coefficient of wavelet decomposition
% Output:
%          Dcoef : the denoising version
N = length(coef);
Dcoef = zeros(N,1);
% The block has the length of three, overlap one
for i = 1:N
    if i == 1
        S = coef(i)^2 + coef(i+1)^2;
    elseif i == N
        S = coef(i-1)^2 + coef(i)^2;
    else
        S = coef(i-1)^2 + coef(i)^2 + coef(i+1)^2;
    end
    if abs(S) < T
        Dcoef(i) = 0;
    else
        Dcoef(i) = coef(i) * (1-T/S);        
    end
end
end
% function y = hardth(x,tau)
% % y:   is the signal after hard-threshold 
% % x:   is the original signal 
% % tau: is the threshold
% y = x.*double((x >= tau)|(-x >= tau));
% end
% 
% function y = softth(x,tau)
% % y:   is the signal after hard-threshold 
% % x:   is the original signal 
% % tau: is the threshold
% y = sign(x).*max(abs(x)-tau,0);
% end


% Made by Zhibin Zhao
% Contact with zhaozhibin@stu.xjtu.edu.cn
% Date: 2017.02.15



