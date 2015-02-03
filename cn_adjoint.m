function P = cn_adjoint(d,PT,u,v,w)
% CN_ADJOINT solve adjoint Bloch equation using adjoint Crank-Nicolson scheme
% M = CN_ADJOINT(D,PT,U,V,W) computes the magnetization vector M starting from
% terminal conditions PT with RF pulse U,V and gradient W. The structure D 
% contains the problem parameters. See
%   C.S. Aigner, C. Clason, A. Rund and R. Stollberger, 
%   Efficient high-resolution RF pulse design applied to simultaneous 
%   multi-slice excitation, 
%   http://math.uni-graz.at/mobis/publications/SFB-Report-2015-001.pdf
%
% February 3, 2015         Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)

steps = length(d.tdis);

I  = eye(3,3);
B1 = d.gamma*d.B1c;

P = zeros(3,d.Nx,steps-1);
parfor z = 1:d.Nx
    B3   = d.gamma*d.G3*d.xdis(z); %#ok<*PFBNS>
    Akp1 = [ -1/d.T2*d.relax,      -w(end)*B3,      -v(end)*B1;...
                   w(end)*B3, -1/d.T2*d.relax,      -u(end)*B1;...
                   v(end)*B1,       u(end)*B1, -1/d.T1*d.relax ]; 
    Pz = zeros(3,steps-1);    
    Pz(:,end) = (I-d.dt/2*Akp1)\PT(:,z,:);
    for k = steps-2:-1:1
        Ak = [ -1/d.T2*d.relax,        -w(k)*B3,        -v(k)*B1;...
                       w(k)*B3, -1/d.T2*d.relax,        -u(k)*B1;...
                       v(k)*B1,         u(k)*B1, -1/d.T1*d.relax ]; 
        Pz(:,k) = (I-d.dt/2*Ak)\((I+d.dt/2*Akp1)*Pz(:,k+1));
        Akp1 = Ak;
    end
    P(:,z,:) = Pz;
end
