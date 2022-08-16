function [DivE,DivH] = div(hx,hy,hz,Nx,Ny,Nz,varargin)

    [DEX,DHX,DEY,DHY,DEZ,DHZ] = fdm.maxwell([hx,hy,hz],[Nx,Ny,Nz],varargin{:});
    
    DivE = [DEX, DEY, DEZ];
    DivH = [DHX, DHY, DHZ];