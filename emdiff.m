 function [varargout] = emdiff(h,N,varargin)
%EMDIFF    Finite difference derivative operators on the Yee grid scheme
%
%   [DEX,DHX] = EMDIFF(h, N) generates finite difference approximation
%   matrices of the first order derivatives for the electric and magnetic
%   fields on the N-D yee grid. The size and dimensions included are taken
%   from the configuration of the h and N arguments as follows,
%
%       [hx],[Nx]                 : generates X derivatives only
%       [hx,hy],[Nx,Ny]           : generates X and Y derivatives
%       [hx,hy,hz],[Nx,Ny,Nz]     : generates X, Y and Z derivatives
%
%   returned are the E and H derivative pairs for each dimension. For
%   example, the following call for a 3-D grid would be,
%
%       [DEX,DHX,DEY,DHY,DEZ,DHZ] = EMDIFF([hx,hy,hz],[Nx,Ny,Nz])
%
%   EMDIFF(__, BC) specifies the numerical boundary conditions to use at 
%   the grid edges. Conditions are specified by coordinate axis, ex. BC = 
%   [xBC,yBC] or for every individual edge, ex. BC=[xminBC,xmaxBC,yminBC,
%   ymaxBC]. The boundary conditions may be one of the following,
%
%       0: Perfect electric conductor (PEC)
%       1: Perfect magnetic conductor (PMC)
%       2: Periodic
%       3: Bloch
%       4: Perfectly matched layers (stretched coordinate PML)
%
%   For Bloch boundaries, an extra argument must be provided to specify the
%   incident wavevector components kin=[kx,ky,kz], where the number of 
%   components must correspond to dimensions.
%
%   For PML boundaries, an extra argument should be provided to configure
%   the number of cells per PML edge. For example, BC=[0,0,3,3] should be
%   followed by PML=[0,0,yminPML,ymaxPML], where yminPML and ymaxPML are
%   the number of cells to use.
%
%   NOTE: to get second order derivative operators from first order ones,
%   simply multiply them together as so (only works in 2D):
%
%       D2EX = DEX * DHX
%       D2EY = DEY * DHY
%       D2HX = -D2EX'
%       D2HY = -D2EX'
%
%%%%%%%%%% Convention, Nx=3, Ny=2 is represented as a 3x2 matrix such that,
%
%       y
%    +---->
%  x | 1 4
%    v 2 5
%      3 6
%
% to plot that matrix, you must transpose it

    oa = 1;
    BC = 0;
    ki = [];
    PML = [];
    
    assert(length(h) == length(N), "h and N must have the same number of elements.")
    
    if nargin > 2
        BC = varargin{1};
    end
    
    assert(length(BC) == 2*length(h), "Not enough boundary conditions provided for each axis.")
    
    if nargin > 3
        if ismember(BC,3)
            ki = varargin{2};
            assert(length(ki) == length(h), "kin must have the same number of components as the number of axes.")
        else
            PML = varargin{2};
            assert(length(PML) == 2*length(h), "PML must have the same number of components as the number of axes.")
        end
    end
    
    if nargin > 4
        PML = varargin{3};
        assert(length(PML) == 2*length(h), "PML must have the same number of components as the number of axes.")
    end

    % fdm weights and offsets
    s = 0:oa; % forward diff
    w = fdm.weights(s);
    
    D = {};

    for i = 1:length(N)
        
        A = sparse(N(i),N(i));
        for c = 1:length(s)
            A = spdiags(w(c) * ones(N(i),1), s(c), A);
        end
        
        switch lower(BC(i))
            case 0 % dirichlet
                % already implemented
                
            case 1
                error("Unimplmeneted feature 'PMC'!")
                
            case 2
                A(end,1) = 1;
                
            case 3
                A(end,1) = exp(1j*ki(i)*N(i)*h(i));
                
            case 4
                % fullfiled below
                
            otherwise
                error("Invalid boundary condition '" + BC + "'")
        end

        for j = 1:i-1
            A = kron(A, eye(N(j)));
        end

        for j = i+1:length(N)
            A = kron(eye(N(j)), A);
        end
        
        D{i} = A;
    end
    
    if ismember(BC,4)
        % apply stretched coordinate PML
        PML(BC < 4) = 0;
        [sx,sy] = pml2d(h*k0,N,PML);
        DEX = fdm.coef(1./sx) * DEX;
        DEY = fdm.coef(1./sy) * DEY;
        DHX = fdm.coef(1./sx) * DHX;
        DHY = fdm.coef(1./sy) * DHY;
    end
    
    for i = 1:numel(D)
        varargout{2*i - 1} =  D{i}  / h(i); % DE
        varargout{2*i - 0} = -D{i}' / h(i); % DH
    end
