function [sx,sy] = scpml2d(Nx,Ny,NPMLx,NPMLy,hx,hy,k0)
% https://web.stanford.edu/group/fan/publication/Shin_JCP_231_3406_2012.pdf

    p = 4;
    
    eta0 = sqrt(4*pi*1e-7 / 8.854187e-12);
    tmax = 16 * (p + 1) / (2*eta0);

    x = zeros(Nx,1);
    x(1:NPMLx) = flip(1:NPMLx) / NPMLx;
    x = flip(x);
    x(1:NPMLx) = flip(1:NPMLx) / NPMLx;
    x = flip(x);
    
    y = zeros(Ny,1);
    y(1:NPMLy) = flip(1:NPMLy) / NPMLy;
    y = flip(y);
    y(1:NPMLy) = flip(1:NPMLy) / NPMLy;
    y = flip(y);
    
    [y,x] = meshgrid(y,x);
    
    tx = tmax/(max(1,NPMLx)*hx) * x.^p;
    ty = tmax/(max(1,NPMLy)*hy) * y.^p;
    
    sx = 1 - 1j*eta0*tx/k0;
    sy = 1 - 1j*eta0*ty/k0;