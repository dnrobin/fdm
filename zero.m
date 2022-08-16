function A = zero(x,varargin)
%ZERO  Identity matrix operator
%
%   I = ZERO(x,...) creates an zero matrix for the N-D grid specified 
%   by the list of coordinate vectors x, y, etc. input.
%
%   I = ZERO([],N) creates a prod(N)xprod(N) zero operator where N may
%   be an array of grid number of points in each axis.

    N = length(x);

    if nargin > 1
        if isempty(x)
            N = varargin{1};
        else
            for i = 1:nargin-1
                N = [N, length(varargin{i})];
            end
        end
    end
    
    A = sparse(prod(N),prod(N));