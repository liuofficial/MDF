function [X] = RLR_MDF(Y, Z, R, BlurD, O, opt)
% A Truncated Matrix Decomposition for Hyperspectral Image Super-Resolution
% Author: Jianjun Liu
% Date: 2019-06-25
O = HSim2mat(O);
max_Y = max(Y(:));
Y = Y ./ max_Y; Z = Z ./ max_Y; O = O ./ max_Y;
% parameters
if ~isfield(opt,'alpha'), alpha = 1; else, alpha = opt.alpha; end
if ~isfield(opt,'lam'), lam = 1e-2; else, lam = opt.lam; end
if ~isfield(opt,'mu'), mu = 0.05; else, mu = opt.mu; end
if ~isfield(opt,'J'), J = 30; else, J = opt.J; end % number of endmembers
if ~isfield(opt,'nC'), nC = 400; else, nC = opt.nC; end 
if ~isfield(opt,'sig'), sig = 25; else, sig = opt.sig; end 
if ~isfield(opt,'maxit'), maxit = 200; else, maxit = opt.maxit; end
if ~isfield(opt,'epsilon'), epsilon = 1e-3; else, epsilon = opt.epsilon; end
if ~isfield(opt,'vc'), vc = 1; else, vc = opt.vc; end
% set
[NW, NH, ~] = size(Z); [~,~,NB] = size(Y);
Yh = zeros(NW,NH,NB);
Yh(1+BlurD.shift:BlurD.ratio:end,1+BlurD.shift:BlurD.ratio:end,:) = Y;
clear Y;
% masking
mask = zeros(NW, NH);
mask(1+BlurD.shift:BlurD.ratio:end,1+BlurD.shift:BlurD.ratio:end) = 1;
mask = repmat(mask, [1, 1, J]);
mask = HSim2mat(mask);
% FFT
FB = cpt_blur_fft(BlurD.K, NW, NH, (size(BlurD.K,1)+1)/2);
FBC = conj(FB);
FBC_BBC = FBC ./ (abs(FB.^2) + 2);
F_BBC = 1 ./ (abs(FB.^2) + 2);
% step 1: initialize endmembers by VCA
Yh = HSim2mat(Yh);
if vc == 1, [A] = cpt_vca(Yh(:, mask(1,:) > 0), J);
else, [A] = cpt_svd(Yh(:, mask(1,:) > 0), J); end
% divide Z into nC superpixels
Z = HSim2mat(Z);
labels = SuperPixelMultiScale(Z, NW, NH, nC, sig);
% initalization
S = zeros(J, NW*NH);
V = cell(3,1);
for i = 1 : 3, V{i} = S; D{i} = 0.*S; end
NU = V;
% step 2: update S
eA = R * A; 
IA = (A'*A + (mu/alpha)*eye(J)) \ eye(J); AtYh = A' * Yh;
IeA = (eA'*eA + mu*eye(J))  \ eye(J); eAtZ = eA' * Z;
for iter = 1 : maxit
    S0 = S;
    S = cpt_conv_fft(V{1}+D{1}, FBC_BBC, NW) + cpt_conv_fft(V{2}+D{2}+V{3}+D{3}, F_BBC, NW);
    NU{1} = cpt_conv_fft(S, FB, NW) - D{1};
    V{1} = IA * (AtYh+(mu/alpha).*NU{1}) .* mask + NU{1}.*(1-mask);
    NU{2} = S - D{2};
    V{2} = IeA * (eAtZ + mu.*NU{2});
    NU{3} = S - D{3};
    V{3} = cpt_sp_svd(NU{3}, lam/mu, labels);
    % update aux
    for i = 1 : 3, D{i} = -NU{i} + V{i}; end
    tol = norm(A*(S0-S),'fro')/norm(A*S0,'fro');
    fprintf('iter: %d, tol, %.4f, %.4f\n', iter, tol, (norm(A*S-O,'fro')^2/numel(O)).^0.5);
    if tol < epsilon, break; end
end
% step 3: update A
eS = HSim2mat( HSBlurDown( HSmat2im(S,NW), BlurD ) );
A = (Yh(:, mask(1,:) > 0)*eS') / (eS*eS'+lam*eye(J));
% step 4: update X
X = A * S;
fprintf('RMSE: %.4f\n', (norm(X-O,'fro')^2/numel(O)).^0.5);
X = HSmat2im(X, NW);
X = X .* max_Y;
end

function [A] = cpt_vca(X, J)
max_vol = 0;
vol = zeros(1, 20);
for idx_VCA = 1:20
    A_aux = VCA(X,'Endmembers',J,'SNR',0,'verbose','off');
    vol(idx_VCA) = abs(det(A_aux'*A_aux));
    if vol(idx_VCA) > max_vol
        A = A_aux;
        max_vol = vol(idx_VCA);
    end
end
end

function [A] = cpt_svd(X, J)
[A,~,~] = svds(X,J);
end

function Y = cpt_sp_svd(X, lam, labels)
Y = zeros(size(X));
label = unique(labels);
nlabel = length(label);
for n = 1 : nlabel
    idx = labels == label(n);
    Xn = X(:,idx);
    [U,S,V] = svd(Xn,'econ');
    S = max(S-lam,0);
    Y(:,idx) = U(:,1:size(S,1))*S*V(:,1:size(S,2))';
end
end

function labels = SuperPixelMultiScale(img, rows, cols, nCs, sigs)
lam = 0.5; dist = 2;
img = normalization(img, 0, 255, 1)';
img = reshape(img, rows, cols, size(img,2));
labels = mex_MSERS(img,nCs,lam,sigs,dist);
end