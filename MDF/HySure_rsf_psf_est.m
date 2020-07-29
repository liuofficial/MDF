function [V, R_est, B_est] = HySure_rsf_psf_est(I_HS, I_MS, ratio, overlap, p, shift)
lambda_R = 1e1;
lambda_B = 1e1;
downsamp_factor = ratio;
ms_bands = size(I_MS,3);
intersection = cell(1,length(ms_bands));
for i = 1:ms_bands
    intersection{i} = overlap(i,:);
end
contiguous = intersection;
% Blur's support: [hsize_h hsize_w]
hsize_h = 2*downsamp_factor-1;
hsize_w = 2*downsamp_factor-1;
blur_center = mod(downsamp_factor+1,2); % to center the blur kernel according to the simluated data
[V, R_est, B_est] = sen_resp_est(I_HS, I_MS, downsamp_factor,...
    intersection, contiguous, p, lambda_R, lambda_B, hsize_h, hsize_w, shift, blur_center);
[nl, nc,~] = size(I_MS);
% Blur kernel
middlel = round((nl+1)/2);
middlec = round((nc+1)/2);
% Blur matrix
B_est = fftshift(B_est);
B_est = B_est(middlel-ratio+1:middlel+ratio-1, middlec-ratio+1:middlec+ratio-1);
end