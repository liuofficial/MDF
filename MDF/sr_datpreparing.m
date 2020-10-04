function [I_REF, I_HS, I_MS, R, BlurD, HS, MS] = sr_datpreparing(I_REF, ratio, SNRh, SNRm, dat_path, spec_rng)
shift = 0; % shift by pixel for down sampling
BlurD.shift = shift; BlurD.ratio = ratio;
[Row, Col, Band] = size(I_REF);
Row = floor(Row/ratio) * ratio; Col = floor(Col/ratio) * ratio;
I_REF = I_REF(1:Row, 1:Col, :);
% Normalize all channels to a given range
[I_REF] = sr_datpruning(I_REF); 
% 5 bands: 1 - pan; 2 - blue; 3 - green; 4 - red; 5 - NIR
ms_bands = 2 : 5; 
load([dat_path 'ikonos_spec_resp.mat']);
[~, valid_ik_bands] = intersect(ikonos_sp(:,1), spec_rng);
no_wa = length(valid_ik_bands);
% Spline interpolation
xx  = linspace(1, no_wa, Band);
x = 1 : no_wa;
R = zeros(5, Band);
for i = 1 : 5 % 1 - pan; 2 - blue; 3 - green; 4 - red; 5 - NIR
    R(i,:) = spline(x, ikonos_sp(valid_ik_bands,i+1), xx);
end
% Use just the predefined bands
R = R(ms_bands,:);
% Spatial degradation (blur)
kersiz = 2*ratio - 1;
sig = (1/(2*(2.7725887)/ratio^2))^0.5;
BlurD.K = fspecial('gaussian',[kersiz kersiz],sig);
I_HS = imfilter(I_REF, BlurD.K, 'circular');
% Add noise
sigmah = sqrt(sum(I_HS(:).^2)/(10^(SNRh/10))/numel(I_HS));
rng(0,'v5uniform');
HS = I_HS;
I_HS = I_HS + sigmah*randn(size(I_HS));
% Down sampling
HS = HS(1+shift:ratio:end,1+shift:ratio:end,:);
I_HS = I_HS(1+shift:ratio:end,1+shift:ratio:end,:);
% Add noise
% sigmah = sqrt(sum(I_HS(:).^2)/(10^(SNRh/10))/numel(I_HS));
% rng(0,'v5uniform');
% I_HS = I_HS + sigmah*randn(size(I_HS));
% Spectral degradation
I_MS = R * reshape(I_REF, [Row*Col Band])';
% Normalize all channels to 1
c = zeros(length(ms_bands),1);
for i=1:length(ms_bands)
    c(i) = max(I_MS(i,:));
    I_MS(i,:) = I_MS(i,:)/c(i);
    R(i,:) =  R(i,:)/c(i);
end
% Add noise
sigmam = sqrt(sum(I_MS(:).^2)/(10^(SNRm/10))/numel(I_MS));
rng(1,'v5uniform');
MS = I_MS;
I_MS = I_MS + sigmam*randn(size(I_MS));
MS = reshape(MS', [Row Col size(MS,1)]);
I_MS = reshape(I_MS', [Row Col size(I_MS,1)]);
end