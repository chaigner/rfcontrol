function [J,G,Xk] = objfun(d,u)
% OBJFUN compute functional value, gradient
% [J,G,XK] = OBJFUN(D,U) computes the value J of the functional 
% to be minimized together with the gradient G in the point U.
% The structure XK contains the necessary information to evaluate the 
% Hessian in U. The structure D contains the problem parameters. See
%   C.S. Aigner, C. Clason, A. Rund and R. Stollberger, 
%   Efficient high-resolution RF pulse design applied to simultaneous 
%   multi-slice excitation, 
%   http://math.uni-graz.at/mobis/publications/SFB-Report-2015-001.pdf
%
% February 3, 2015         Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)

% zero padding of control to readout time
u = [u; zeros(d.Nt-1-d.Nu,1)];

% solve state equation
M = cn_bloch(d,d.M0,u,d.v,d.w);

% residual
res = M(:,:,end) - d.Md;

% objective function
J = 1/2*d.dx*norm(res(:))^2 + d.alpha/2*d.dt*norm(u(1:d.Nu))^2;

% gradient
if nargout>1
    P = cn_adjoint(d,res,u,d.v,d.w);
    N = 1/2*(M(:,:,1:end-1) + M(:,:,2:end));
    G = d.alpha*u(1:d.Nu) + d.gamma*d.B1c*...
          d.dx*sum(squeeze(N(3,:,1:d.Nu).*P(2,:,1:d.Nu) - N(2,:,1:d.Nu).*P(3,:,1:d.Nu)))';
end

% Hessian information
if nargout>2
    Xk.N = N;
    Xk.P = P;
    Xk.u = u;
end
