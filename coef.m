function C = coef(c)
%COEF  Constant coefficient matrix operator
%
%   C = COEF(c) creates an NxN diagonal matrix from the coefficent array
%   placing the N coefficients on the main diagonal.

    C = diag(sparse(c(:)), 0);