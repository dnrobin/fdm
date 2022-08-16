function s = stencil(n, varargin)
%STENCIL    Nearest grid neighbour offsets to calculate finite differences
%
%   s = STENCIL(n) returns an array of nearest neighbour offsets of a grid
%   cell generated to compute a central finite difference approximation to
%   the n'th order derivative with the minimal order of accuracy needed.
%
%   s = STENCIL(__, acc) specifies the order of accuracy (default 1) after 
%   the previous syntax.
%
%   s = STENCIL(__, scheme) optionally specify the finite difference scheme
%   as one of 'forward', 'backward', (default) 'central' after any of the
%   previous syntaxes. Note that for central fdm, the order of accuracy is
%   always rounded to the nearest even number, ie., acc = 2, 4, 6, 8, etc.

    acc = 1;
    scheme = 'central';
    
    if nargin > 1
        if isnumeric(varargin{1})
            acc = varargin{1};
            varargin(1) = [];
        end
        if ~isempty(varargin)
            scheme = varargin{1};
        end
    end
    
    assert(acc > 0, "Order of accuracy must be greater than zero")

    m = n + acc - 1;
    
    switch lower(scheme)
        case 'forward'
            s = 0:m;
        case 'backward'
            s = -m:0;
        case 'central'
            d = max(1, floor((m + 1) / 2));
            s = -d:d;
        otherwise
            error("Invalid scheme argument '" + scheme + "'")
    end
