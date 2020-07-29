function [A] = cpt_conv_fft(X, FK, NW)
[NB, N] = size(X);
NH = N/NW;
A = reshape(real(ifft2(fft2(reshape(X', NW, NH, NB)).*repmat(FK,[1, 1, NB]))), NW*NH, NB)';
end