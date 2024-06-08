function [S0,b0,c0,txt] = get_obj(theta,d,k,xi,M0,h,einc_vect,obj_no,x,y,Nx,Ny,Mx,My,einc_pw_fun,const,bg_idx_r,varargin)

dxdy = h^2;
KK = size(M0,1);

switch obj_no
    case 1 % objective matrices (Pabs)
        S0 = (1/2/k^3)*Im(xi*(M0'*M0))*dxdy;
        b0 = zeros(KK,1);
        c0 = 0;
        txt = 'Pabs';
    case 2 % objective matrices (Pscat)
        S0 = (1/2/k)*Im(M0)*dxdy;
        b0 = zeros(KK,1);
        c0 = 0;
        txt = 'Pscat';
    case 3 % objective matrices (Pext)
        S0 = zeros(size(M0));
        b0 = -1j*dxdy/2/k*M0'*einc_vect;
        c0 = 0;
        txt = 'Pext';
    case 4 % objective matrices (|Es|^2)
        S0 = M0'*M0*dxdy;
        b0 = zeros(KK,1);
        c0 = 0;
        txt = '|Es|^2';
        
    case 5
        [G0,einc0] = focusing(k,x,y,Nx,Ny,Mx,My,einc_pw_fun,const);
        
        G0 = -k^2*dxdy*G0; % row vector
        G0(bg_idx_r) = 0;
        
        S0 = (1/k^4)*(M0'*G0')*(G0*M0)*dxdy;
        b0 = (-1/k^2)*2*M0'*G0'*einc0*dxdy;
        c0 = abs(einc0)^2*dxdy;
        txt = '|E(target)|^2';             
        
    case 6
        
        [G0,einc0] = focusing(d,k,x,y,Nx,Ny,Mx,My,einc_pw_fun,const);
        
        G0 = -k^2*dxdy*G0; % row vector
        G0(bg_idx_r) = 0;
        
        S0 = zeros(size(M0));
        b0 = (-1/k^2)*exp(-1i*theta)*M0'*G0';
        c0 = real(einc0*exp(1i*theta));
        txt = '|E(target)|^2 theta'; 
end

if ~isempty(varargin) && strcmpi('sparse',varargin{1})
    S0 = sparse(S0);
    b0 = sparse(b0);
    c0 = sparse(c0);    
end

end