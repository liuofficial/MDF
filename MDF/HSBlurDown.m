function X = HSBlurDown(X, BlurD)
X = imfilter(X, BlurD.K, 'circular');
X = X(1+BlurD.shift:BlurD.ratio:end,1+BlurD.shift:BlurD.ratio:end,:);
end