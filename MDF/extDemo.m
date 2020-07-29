% For more information, please see "D. Shen, J. Liu, Z. Xiao, J. Yang and L. Xiao, "A Twice Optimizing Net With Matrix Decomposition for Hyperspectral and Multispectral Image Fusion," in IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 13, pp. 4095-4110, 2020"
% files: data/info.mat (R, B); data/CAVE/21.mat...32.mat (HS, LRHS, HRMS)
clear; clc;
data_path = 'data/';
vca_fea = 10; svd_fea = 31;
ratio = 8;
J = 10; % number of endmembers
mat = load([data_path 'info.mat']);
R_true = mat.R; B_true = mat.B;
qtys = zeros(6, 12, 6);
BlurD.shift = 0;
BlurD.ratio = ratio;
for i = 1 : 12
    disp(['=======================' num2str(i) '==========================='])
    mat = load([data_path 'CAVE/' num2str(i+20) '.mat']);
    I_REF = double(mat.HS); I_HS = double(mat.LRHS); I_MS = double(mat.HRMS);
    % estimate R and B
    p = J;
    overlap = repmat(1:size(I_HS,3),[size(I_MS,3) 1]);
    [V, R_est, B_est] = HySure_rsf_psf_est(I_HS, I_MS, ratio, overlap, p, BlurD.shift);
    disp('---CAVE, Noise Input Data: B/R est---');
    R = R_est; BlurD.K = B_est;
    disp('---VCA: LSMDF_ver---');
    opt = struct;
    opt.J = vca_fea;
    opt.vc = 1;
    Out = LSMDF2(I_HS, I_MS, R, BlurD, opt);
    [~,qty] = HSQualityIndices(Out,I_REF,ratio);
    qtys(:,i,1) = qty;
    disp('---SVD: LSMDF_ver---');
    opt = struct;
    opt.J = svd_fea;
    opt.vc = 0;
    Out = LSMDF2(I_HS, I_MS, R, BlurD, opt);
    [~,qty] = HSQualityIndices(Out,I_REF,ratio);
    qtys(:,i,2) = qty;
    
    disp('---CAVE, Noise Input Data: B/R true---');
    R = R_true; BlurD.K = B_true;
    disp('---VCA: LS_MDF_ver---');
    opt = struct;
    opt.J = vca_fea;
    opt.vc = 1;
    Out = LSMDF2(I_HS, I_MS, R, BlurD, opt);
    [~,qty] = HSQualityIndices(Out,I_REF,ratio);
    qtys(:,i,3) = qty;
    disp('---SVD: LSMDF_ver---');
    opt = struct;
    opt.J = svd_fea;
    opt.vc = 0;
    Out = LSMDF_ver(I_HS, I_MS, R, BlurD, opt);
    [~,qty] = HSQualityIndices(Out,I_REF,ratio);
    qtys(:,i,4) = qty;
    
    disp('---CAVE, True Input Data: B/R true---');
    R = R_true; BlurD.K = B_true;
    HS = HSBlurDown(I_REF, BlurD);
    MS = R * HSim2mat(I_REF);
    MS = HSmat2im(MS, size(I_REF,1));
    disp('---VCA: LSMDF_ver---');
    opt = struct;
    opt.J = vca_fea;
    opt.vc = 1;
    Out = LSMDF2(HS, MS, R, BlurD, opt);
    [~,qty] = HSQualityIndices(Out,I_REF,ratio);
    qtys(:,i,5) = qty;
    disp('---SVD: LSMDF_ver---');
    opt = struct;
    opt.J = svd_fea;
    opt.vc = 0;
    Out = LSMDF2(HS, MS, R, BlurD, opt);
    [~,qty] = HSQualityIndices(Out,I_REF,ratio);
    qtys(:,i,6) = qty;
end