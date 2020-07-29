function [out, outall] = HSPSNR(ref,tar,mask)
%--------------------------------------------------------------------------
% Peak signal to noise ratio (PSNR)
%
% USAGE
%   out = PSNR(ref,tar,mask)
%
% INPUT
%   ref : reference HS data (rows,cols,bands)
%   tar : target HS data (rows,cols,bands)
%   mask: binary mask (rows,cols) (optional)
%
% OUTPUT
%   out : PSNR (scalar)
%
%--------------------------------------------------------------------------
[~,~,bands] = size(ref);
if nargin == 2
    ref = reshape(ref,[],bands);
    tar = reshape(tar,[],bands);
    msr = mean((ref-tar).^2,1);
    max2 = max(ref,[],1).^2;
    
    psnrall = 10*log10(max2./msr);
    outall.all = psnrall;
    outall.ave = mean(psnrall);
else
    ref = reshape(ref,[],bands);
    tar = reshape(tar,[],bands);
    msr = mean((ref(mask~=0,:)-tar(mask~=0,:)).^2,1);
    max2 = max(ref,[],1).^2;
    
    psnrall = 10*log10(max2./msr);
    outall.all = psnrall;
    outall.ave = mean(psnrall);
end
out = outall.ave;
end