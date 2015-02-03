function u = tr_newton(d,param,objfun,applyHess,u0)
% TR_NEWTON trust-region CG-Newton method
% U = TR_NEWTON(D,PARAM,OBJFUN,APPLYHESS,U0) computes the optimal control U
% using a trust-region CG-Newton method. The functional to be minimized is
% specified using the function handles OBJFUN, which evaluates functional
% and gradient, and APPLYHESS, which computes the action of the Hessian on
% a given direction. The structure PARAM contains the necessary parameters
% for the trust-region Newton method, while the structure D contains the 
% problem parameters. U0 is the initial guess for the control U. See
%   C.S. Aigner, C. Clason, A. Rund and R. Stollberger, 
%   Efficient high-resolution RF pulse design applied to simultaneous 
%   multi-slice excitation, 
%   http://math.uni-graz.at/mobis/publications/SFB-Report-2015-001.pdf
%
% February 3, 2015         Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)

maxit  = param.maxit;       % maximum number of TR iteration (20)
reltol = param.reltol;      % terminate if gradient reduced by (1e-4)
abstol = param.abstol;      % terminate if gradient smaller (1.2e-7)
rho    = param.rho;         % initial trust region radius (1)
maxrad = param.maxrad;      % maximal trust region radius (70)
sig1   = param.sig1;        % increase trad if dJa/dJm > sig1 (0.7)
sig2   = param.sig2;        % no change if sig2 < dJa/dJm < sig2 (0.25)
sig3   = param.sig3;        % increase trad if dJa/dJm > sig3 (0.7)
q      = param.q;           % factor for radius change (4)
ip     = param.ip;          % inner product for CG
cgtol  = param.cgtol;       % desired reduction of residual in CG
cgits  = param.cgits;       % maximum number of CG iterations

[J,g,Xk] = objfun(d,u0);
nrG0 = sqrt(ip(g,g));   it = 0;    u = u0;

fprintf('it \tJ \t\t|g| \t\tflag \trho \t\tdJa/dJm \tcgits\n')
fprintf([num2str(it)  '\t' num2str(J,'%1.3e') '\t' ...
    num2str(nrG0,'%1.3e') ' \n']);

for it = 1:maxit
    % minimize quadratic model
    Hmult = @(du) applyHess(d,Xk,du);
     
    [du,flag,cgit] = tr_cg(Hmult,g,ip,rho,cgtol,cgits);
    
    % test if control is updated
    dJa = J - objfun(d,u+du);                   % actual reduction in J
    dJm = -(1/2*ip(du,Hmult(du)) + ip(du,g));   % predicted reduction in J
    Jratio = dJa/dJm;               % ratio of real and predicted decrease
    if dJa > eps && Jratio > sig1   % actual reduction and model good
        u = u + du;                 % accept step
        [J,g,Xk] = objfun(d,u);
    else
        fprintf('R');
    end
    
    % test if radius is updated
    if dJa > eps && abs(Jratio-1) <= 1-sig3   % step accepted, model good
        rho = min(q*rho,maxrad);              % increase radius
    elseif dJa <= eps                         % step rejected, no decrease
        rho = 1/q*rho;                        % decrease radius
    elseif  Jratio < sig2                     % model bad
        rho = 1/q*rho;                        % decrease radius
    end
   
    nrG = sqrt(ip(g,g));
    fprintf([num2str(it) '\t' num2str(J,'%1.3e') '\t' ...
        num2str(nrG,'%1.3e') '\t' num2str(flag) '\t'  ...
        num2str(rho,'%1.3e') '\t' num2str(Jratio,'%1.3e') '\t' ...
        num2str(cgit) ' \n']);
    
    if (nrG < reltol*nrG0 || nrG < abstol )      % tolerance reached
        break
    end
      
end
