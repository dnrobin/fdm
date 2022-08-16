function w = weights(x, varargin)
%WEIGHTS    Weights coefficients to calculate finite differences
%
%   w = WEIGHTS(x) returns an array of weights for the given offset grid
%   coordinates to calculate the finite difference approximation to a first
%   order derivative on those grid points.
%
%   w = WEIGHTS(__, n) specifies the derivative order (default 1).

    n = 1;
    N = length(x);
    
    if nargin > 1
        n = varargin{1};
    end

    if n + 1 > N
        error('Number of grid points must be at least n + 1')
    end

    S = zeros(N);
    for i = 1:N
        S(:,i) = x(:).^(i-1);
    end
    
    b = zeros(N,1);
    b(n + 1, :) = factorial(n);
    
    w = S' \ b(:);
end