function [du,flag,it] = tr_cg(H,g,ip,trad,tol,maxit)
% TR_CG trust-region conjugate gradient iteration
% [DU,FLAG,IT] = TR_CG(H,G,IP,TRAD) solves the Newton step HDU = -G using
% Steihaug's trust-region conjugate gradient method, where H is a function
% handle for the evaluation of H(X), IP is a function handle for computing
% the inner product used in the CG iteration, and TRAD is the radius of the
% trust region. The convergence FLAG returns:
%    0 TRCG converged to the desired tolerance TOL within IT iterations
%    1 TRCG iterated MAXIT times but did not converge
%    2 TRCG terminated because the iterate left the trust region
%    3 TRCG terminated because negative curvature was encountered
% See
%   C.S. Aigner, C. Clason, A. Rund and R. Stollberger, 
%   Efficient high-resolution RF pulse design applied to simultaneous 
%   multi-slice excitation, 
%   http://math.uni-graz.at/mobis/publications/SFB-Report-2015-001.pdf
%
% February 3, 2015         Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)

du = 0*g;
r  = -g;
p  = r;
nr = ip(r,r);    nr0 = sqrt(nr);    it = 1;
while true
    Hp  = H(p);
    pHp = ip(p,Hp);
    if pHp < eps                        % negative curvature
        tau = dist2bdy(du,p,trad,ip);   % go to boundary
        du  = du + tau*p;
        flag = 3; break;
    end
    
    al = nr/pHp;
    if ip(du+al*p,du+al*p) >= trad^2    % step too large
        tau = dist2bdy(du,p,trad,ip);   % go to boundary
        du = du + tau*p;
        flag = 2; break;
    end
    
    du  = du + al*p;
    r   = r - al*Hp;
    nrk = ip(r,r);
    
    if nrk < tol*nr0^(1.3)              % norm of residual small enough
        flag = 0; break;
    elseif it == maxit                  % too many iterations
        flag = 1; break;
    end
    
    p  = r + (nrk/nr)*p;
    nr = nrk;
    it = it+1;
end


function tau = dist2bdy(du, p, trad, ip)
% find distance to trust-region boundary from du in direction p
dd = ip(p,p);
xd = ip(du,p);
xx = ip(du,du);
ss = trad*trad;

det = xd*xd + dd * (ss - xx);
tau = (ss - xx) / (xd + sqrt(det));
