function [X, A, S] = LSMDF(Y, Z, R, BlurD, opt)
max_Y = max(Y(:));
Y = Y ./ max_Y; Z = Z ./ max_Y;
if ~isfield(opt,'lam'), lam = 1e-2; else, lam = opt.lam; end
if ~isfield(opt,'J'), J = 30; else, J = opt.J; end
if ~isfield(opt,'vc'), vc = 1; else, vc = opt.vc; end
[NW, ~, ~] = size(Z);
Y = HSim2mat(Y); Z = HSim2mat(Z);
if vc == 1, [A] = cpt_vca(Y, J);
else, [A] = cpt_svd(Y, J); end
eA = R * A;
S = (eA'*eA+lam*eye(J)) \ (eA'*Z);
eS = HSim2mat( HSBlurDown( HSmat2im(S,NW), BlurD ) );
A = (Y*eS') / (eS*eS'+lam*eye(J));
X = A * S;
X = HSmat2im(X,NW);
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