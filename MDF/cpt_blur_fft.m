function FB = cpt_blur_fft(B, NW, NH, ratio)
halfW = round((NW+1)/2);
halfH = round((NH+1)/2);
FB = zeros(NW,NH);
FB(halfW-ratio+1:halfW+ratio-1, halfH-ratio+1:halfH+ratio-1) = B;
FB = ifftshift(FB);
FB = fft2(FB);
end