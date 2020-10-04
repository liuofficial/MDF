function [data, maxV, minV] = sr_datpruning(img)
rdown = 0.025; rup = 0.975;
% rdown = 1e-4; rup = 1 - 1e-3;
[data, maxV, minV] = DPTailor(img, rdown, rup);
end

function [img, lmax, lmin] = DPTailor(dat, rdown, rup)
img = dat;
[lmax, lmin] = linear(dat, rdown, rup);
img(dat<lmin) = lmin; img(dat>lmax) = lmax;
img = img / lmax;
end

function [lmax, lmin] = linear(X, a, b)
X = X(:);
x = sort(X);
L = length(x);
lmin = x(ceil(L*a));
lmax = x(ceil(L*b));
end