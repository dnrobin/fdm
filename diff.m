function [varargout] = diff(x, varargin)
% DIFF  Finite difference approximation to derivatives
%
%   D = DIFF(x) returns a first order, one-dimensional finite difference
%   derivative operator matrix for the grid coordinates specified by x.
%
%   [DX,DY,...] = DIFF(x,y,...) returns first order derivative operators
%   on an N-D grid of points given by x, y, ... coordinate vectors. The 
%   grid size does not need to be square.
%
%   [DX,...] = DIFF([], h, N) same as previous but on a uniform grid with 
%   spacing h and number of points N. For N-D grids, the inputs should
%   be arrays such that, h=[hx,hy,...], N=[Nx,Ny,...]. The number of
%   returned arguments will match ne number of input dimensions.
%
%   DIFF(__, n) optionally specify the derivative order n after any of the 
%   previous syntaxes (default n = 1).
%
%   DIFF(__, NAME, VALUE) the function also takes the following name/value
%   pair arguments,
%       'Accuracy' - order of accuracy (default 2)
%       'Stencil'  - custom stencil points array (see fdm.stencil)
%       'Scheme'   - one of 'forward', 'backward', (default) 'central'
%       'Edge'     - improve edge accuracy, 'auto' or 'none'

    accuracy = 1;
    stencil  = [];
    scheme   = 'central';
    edge     = 'none';
    
    n = 1;
    h = [];
    N = length(x);
    
    if nargin > 1
        if isempty(x)
            
            assert(nargin > 2, "Expecting grid number of points N along with spacing argument")
            h = varargin{1};
            N = varargin{2};
            varargin(1:2) = [];
            
            assert(length(h) == length(N), "Size of array h must match that of N")
        else
            
            error("Unimplemented feature: custum grid")
            
            % TODO: check derivative order n after grid input
        end
    else
        error("Unimplemented feature: custum grid")
    end

    if ~isempty(varargin)
        if isnumeric(varargin{1})
            n = varargin{1};
            varargin(1) = [];
        end
    end
    
    for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'accuracy'
                accuracy = varargin{i + 1};
            case 'stencil'
                stencil = varargin{i + 1};
            case 'scheme'
                scheme = varargin{i + 1};
            case 'edge'
                edge = varargin{i + 1};
            otherwise
                error("Invalid name/value pair argument '" + varargin{i} + "'")
        end
    end
    
    % generate stencil
    if ~isempty(stencil)
        s = stencil(:)';
    else
        s = fdm.stencil(n, accuracy, scheme);
    end

    if numel(s) > N(1)
%         error("Too few grid points for the required order of accuracy")
    end
    
    % fdm weights
    w = fdm.weights(s, n);

    % build the kernel
    A = spalloc(N(1),N(1),N(1)*length(s));
    for i = 1:length(s)
        A = spdiags(w(i) * ones(N(1),1), s(i), A);
    end

    % improve accuracy at grid edges
    if edge == "auto"
        s = fdm.stencil(n, accuracy, 'forward');
        for i = 1:length(s)
            c = s - s(i);
            b = c + i;
            if s(i) > 0
                b = b + (N(1) - length(c));
            end
            A(b(i), b) = fdm.weights(c, n);
        end
    end

    % build derivative operators
    D = {A};

    for j = 2:length(N)

        B = {};

        for i = 1:length(D)
            B{end + 1} = kron(speye(N(j)), D{i});
        end
        B{end + 1} = kron(D{i}, speye(N(j)));

        D = B;
    end

    % fit to grid spacing
    for i = 1:numel(D)
        varargout{i} = D{i} / h(i)^n;
    end

    