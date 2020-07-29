function [Out, val] = HSQualityIndices(I_HS,I_REF,ratio,border)
[rows,cols,bands] = size(I_REF);

if nargin < 4
    border = 0;
end
if border == 1
    %Remove border from the analysis
    I_HS  = I_HS(ratio+1:rows-ratio,ratio+1:cols-ratio,:);
    I_REF = I_REF(ratio+1:rows-ratio,ratio+1:cols-ratio,:);
end

Out.rmse = HSRMSE(I_REF,I_HS);
Out.psnr = HSPSNR(I_REF,I_HS);
Out.sam = HSSAM(I_REF,I_HS);
Out.cc = HSCC(I_REF,I_HS);
Out.ergas = HSERGAS(I_REF,I_HS,ratio);
Out.uiqi = HSUIQI(I_REF,I_HS);
%Out.q2n = HSq2n(I_REF, I_HS, 32, 32);
%Out.ssim = HSSSIM(I_REF, I_HS, 0, 0);
%Out.dd = HSDD(I_REF, I_HS);

disp(['RMSE(0) : ' num2str(Out.rmse)]);
disp(['PSNR(inf) : ' num2str(Out.psnr)]);
disp(['SAM(0)  : ' num2str(Out.sam)]);
disp(['CC(1)   : ' num2str(Out.cc)]);
disp(['ERGAS(0): ' num2str(Out.ergas)]);
disp(['UIQI(1)  : ' num2str(Out.uiqi)]);
%disp(['Q2n(1)  : ' num2str(Out.q2n)]);

val = [Out.rmse;Out.psnr;Out.sam;Out.cc;Out.ergas;Out.uiqi];
end
