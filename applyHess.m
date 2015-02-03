function Hdu = applyHess(d,Xk,du)
% APPLYHESS compute application of Hessian
% HDU = APPLYHESS(D,XK,DU) computes the action HDU of the Hessian in 
% direction DU. The structure XK contains the point in which the Hessian
% is evaluated. The structure D contains the problem parameters. See
%   C.S. Aigner, C. Clason, A. Rund and R. Stollberger, 
%   Efficient high-resolution RF pulse design applied to simultaneous 
%   multi-slice excitation, 
%   http://math.uni-graz.at/mobis/publications/SFB-Report-2015-001.pdf
%
% February 3, 2015         Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)


N = Xk.N;     P = Xk.P;    u = Xk.u;  v = d.v;    w = d.w;

steps = length(d.tdis);
I     = eye(3,3);
B1    = d.gamma*d.B1c;

Hdu = d.alpha*du;

% zero padding to readout time
du = [du;zeros(d.Nt-1-d.Nu,1)];

parfor z = 1:d.Nx
    Nz = N(:,z,:);    Pz = P(:,z,:);
    B3 = d.gamma*d.G3*d.xdis(z); %#ok<*PFBNS>
    % solve linearized state equation
    dMz  = zeros(3,steps);
    for k = 2:steps
        bk = [0; B1*Nz(3,k-1)*du(k-1); -B1*Nz(2,k-1)*du(k-1)];
        Ak = [ -1/d.T2*d.relax,       w(k-1)*B3,       v(k-1)*B1;...
                    -w(k-1)*B3, -1/d.T2*d.relax,       u(k-1)*B1;...
                    -v(k-1)*B1,      -u(k-1)*B1, -1/d.T1*d.relax ];
        dMz(:,k) = (I-d.dt/2*Ak) \ ((I+d.dt/2*Ak)*dMz(:,k-1) + d.dt*bk);
    end
    dNz = (dMz(:,1:d.Nu)+dMz(:,2:d.Nu+1))/2;
    dq = dMz(:,end);
    
    % solve linearized adjoint equation
    Akp1 = Ak';
    bkp1 = [0; -B1*Pz(3,end)*du(end); B1*Pz(2,end)*du(end)];
    dPz  = zeros(3,steps-1);   
    dPz(:,end) = (I-d.dt/2*Akp1)\(dq+d.dt/2*bkp1);
    for k = steps-2:-1:1
        bk = [0; -B1*Pz(3,k)*du(k); B1*Pz(2,k)*du(k)];
        Ak = [ -1/d.T2*d.relax,        -w(k)*B3,        -v(k)*B1;...
                       w(k)*B3, -1/d.T2*d.relax,        -u(k)*B1;...
                       v(k)*B1,         u(k)*B1, -1/d.T1*d.relax ]; 
        dPz(:,k) = (I-d.dt/2*Ak) \ ...
            ((I+d.dt/2*Akp1)*dPz(:,k+1) + d.dt/2*(bk+bkp1));
        Akp1 = Ak;    bkp1 = bk;
    end
    dPz = dPz(:,1:d.Nu);
    
    % action of Hessian
    Hdu = Hdu + B1*d.dx*(dNz(3,:).* Pz(2,1:d.Nu) - dNz(2,:).* Pz(3,1:d.Nu) + ...
                          Nz(3,1:d.Nu).*dPz(2,:) -  Nz(2,1:d.Nu).*dPz(3,:))';
end
