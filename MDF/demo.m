%% A Truncated Matrix Decomposition for Hyperspectral Image Super-Resolution, TIP, 2020
% Note: B(bands) >> J(endmembers)
% For the HSI such as CAVE, please see extDemo.m
clear; clc;
ratio = 5;
SNRh = 30; % SNR (in dB) for HSI
SNRm = 40; % SNR (in dB) for MSI
J = 30; % number of endmembers
load(['paviaU.mat']); spec_rng = 430 : 860; I_REF = I_REF(1:580,:,:);
[I_REF, I_HS, I_MS, R, BlurD] = sr_datpreparing(I_REF, ratio, SNRh, SNRm, '', spec_rng);
%% estimate R and BlurD
overlap = repmat(1:size(I_HS,3),[size(I_MS,3) 1]);
p = J;
[V, Rest, Best] = HySure_rsf_psf_est(I_HS, I_MS, ratio, overlap, p, BlurD.shift);
%% R and BlurD
R = Rest; BlurD.K = Best;
%% paras optimization
opt = struct;
opt.J = J;
opt_rlr = struct;
opt_rlr.nC = 200; opt_rlr.sig = 25; opt_rlr.J = J;
%% algorithms
t0 = tic;
Out = LSMDF(I_HS, I_MS, R, BlurD, opt);
disp(toc(t0));
[~,qty] = HSQualityIndices(Out,I_REF,ratio);
t0 = tic;
Out = RLR_MDF(I_HS, I_MS, R, BlurD, I_REF, opt_rlr);
disp(toc(t0));
[~,qty] = HSQualityIndices(Out,I_REF,ratio);

