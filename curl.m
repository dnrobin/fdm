function [CurlE,CurlH] = curl(hx,hy,hz,Nx,Ny,Nz,varargin)
    [DEX,DHX,DEY,DHY,DEZ,DHZ] = fdm.maxwell([hx,hy,hz],[Nx,Ny,Nz],varargin{:});
    
    Z = zeros(Nx*Ny*Nz);
    
    Ce = [
        Z, -DEZ, DEY
        DEZ, Z, -DEX
        -DEY, DEZ, Z
    ];

    Ch = [
        Z, -DHZ, DHY
        DHZ, Z, -DHX
        -DHY, DHZ, Z
    ];