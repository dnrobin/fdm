function [sx,sy] = upml2d(NGRID,NPML,a,t,p)
% CALCPML2D Calculate the PML parameters on a 2D grid
%
% [sx,sy] = CALCPML2D(NGRID,NPML);
%
% This MATLAB function calculates the PML parameters sx and sy
% to absorb outgoing waves on a 2D grid.
%
% Input Arguments
% =================
% NGRID Array containing the number of points in the grid
% = [ Nx Ny ]
% NPML Array containing the size of the PML at each boundary
% = [ Nxlo Nxhi Nylo Nyhi ]
%
% Output Arguments
% =================
% sx,sy 2D arrays containing the PML parameters on a 2D grid
    
    Nx = NGRID(1);
    Ny = NGRID(2);
    
    Nxlo = NPML(1);
    Nxhi = NPML(2);
    Nylo = NPML(3);
    Nyhi = NPML(4);

    eta0 = sqrt(4*pi*1e-7 / 8.854187e-12);
    
    if nargin < 3
        a = 3;
    end
    if nargin < 4
        t = 1;
    end
    if nargin < 5
        p = 3;
    end

    x = zeros(Nx,1);
    x(1:Nxlo) = flip(1:Nxlo) / Nxlo;
    x = flip(x);
    x(1:Nxhi) = flip(1:Nxhi) / Nxhi;
    x = flip(x);
    
    y = zeros(1,Ny);
    y(1:Nylo) = flip(1:Nylo) / Nylo;
    y = flip(y);
    y(1:Nyhi) = flip(1:Nyhi) / Nyhi;
    y = flip(y);
    
    [y,x] = meshgrid(y,x);
    
    ax = 1 + a * x.^p;
    ay = 1 + a * y.^p;
    
    tx = t * sin(pi/2 * x).^2;
    ty = t * sin(pi/2 * y).^2;
    
    sx = ax .* (1 + 1j*eta0*tx);
    sy = ay .* (1 + 1j*eta0*ty);
