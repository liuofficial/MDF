function [A] = HSmat2im(X, row)
%mat2im - converts a matrix to a 3D image
[p, n] = size(X);
nc = n/row;
A = reshape(X', row, nc, p);
end