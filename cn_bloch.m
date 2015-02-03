function M = cn_bloch(d,M0,u,v,w)
% CN_BLOCH solve Bloch equation using Crank-Nicolson scheme
% M = CN_BLOCH(D,M0,U,V,W) computes the magnetization vector M starting
% from initial conditions M0 with RF pulse U,V and gradient W. The
% structure D contains the problem parameters. See
%   C.S. Aigner, C. Clason, A. Rund and R. Stollberger, 
%   Efficient high-resolution RF pulse design applied to simultaneous 
%   multi-slice excitation, 
%   http://math.uni-graz.at/mobis/publications/SFB-Report-2015-001.pdf
%
% February 3, 2015         Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)

steps = length(d.tdis);

b  =[0;0;d.M0c]/d.T1*d.relax;
I  = eye(3,3);
B1 = d.gamma*d.B1c;

M = zeros(3,d.Nx,steps);
parfor z = 1:d.Nx
    B3 = d.gamma*d.G3*d.xdis(z); %#ok<*PFBNS>
    Mz = zeros(3,steps); Mz(:,1) = M0(:,z);
    for k = 2:steps
        Ak = [ -1/d.T2*d.relax,       w(k-1)*B3,       v(k-1)*B1;...
                    -w(k-1)*B3, -1/d.T2*d.relax,       u(k-1)*B1;...
                    -v(k-1)*B1,      -u(k-1)*B1, -1/d.T1*d.relax ];
        Mz(:,k) = (I-d.dt/2*Ak)\((I+d.dt/2*Ak)*Mz(:,k-1) +  d.dt*b);
    end
    M(:,z,:) = Mz;
end
