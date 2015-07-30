%TEST_MULTI test script for multi-slice excitation
% This m-file computes the optimal RF pulse for the simultaneous multi-
% slice problem using the approach and the parameters described in the paper
%   C.S. Aigner, C. Clason, A. Rund and R. Stollberger, 
%   Efficient high-resolution RF pulse design applied to simultaneous 
%   multi-slice excitation, 
%   http://math.uni-graz.at/mobis/publications/SFB-Report-2015-001.pdf
%
% February 3, 2015,  V1.0 original version 
%    July 30, 2015,  V1.1 changed to a default flip angle of 90 degrees 
%
%                          Christoph S. Aigner (christoph.aigner@tugraz.at)
%                          Christian Clason (christian.clason@uni-due.de)
%                          Armin Rund (armin.rund@uni-graz.at)


%% set parameters
% space discretization
d.a    = 0.5;                     % domain border in m
d.z    = 0.0025;                  % half slice thickness in m
d.Nx   = 5001;                    % total number of spatial points
d.xdis = linspace(-d.a,d.a,d.Nx); % spatial running variable
d.dx   = d.xdis(2)-d.xdis(1);     % spatial grid size

% time discretization
load('z_grad_thk2_dt5.mat');     % slice selective gradient shape
d.T     = 13.92;                 % optimization time in ms 
d.Nt    = size(z_grad,1)+1;      % total number of temporal points
d.tdis  = linspace(0,d.T,d.Nt);  % temporal running variable
d.dt    = d.tdis(2) - d.tdis(1); % temporal grid size
d.Nu    = 512;                   % number of temporal control points

% model parameters
d.gamma = 267.51;       % gyromagnetic ratio 42.57*2*pi in MHz/T
d.T1    = 102;          % longitudinal relaxation time in ms
d.T2    = 81;           % transversal relaxation time in ms
d.B0    = 3000;         % static magnetic field strength in mT
d.M0c   = 1;            % normalized equilibrium magnetization
d.B1c   = 0.5e-2;       % weighting for the RF amplitude [u*1e3*d.B1c] = muT
d.G3    = 0.25;         % weighting for the z-Gradient in mT
d.relax = 0;            % 0=without relaxation, 1=with relaxation
d.u0 = zeros(d.Nu,1);   % RF initial guess
d.v  = zeros(d.Nt-1,1); % fixed with zeros
d.w  = z_grad;          % fixed with external shape
d.alpha = 1e-4;         % control costs for u (SAR)

% initial magnetization
d.M0    = d.M0c*repmat([0;0;1],1,d.Nx); 

% TR-CG-Newton parameters
tr.maxit  = 5;        % maximum number of TR Newton iterations
tr.reltol = 1e-4;      % relative tolerance for gradient norm in Newton
tr.abstol = 1.2e-7;    % absolute tolerance for gradient norm in Newton
tr.rho    = 1;         % initial trust region radius
tr.maxrad = 2;         % maximal trust region radius
tr.sig1   = 0.03;      % decrease trad if dJa/dJm < sig1        
tr.sig2   = 0.25;      % do not change if sig1 < dJa/dJm < sig2
tr.sig3   = 0.7;       % increase trad if dJa/dJm > sig3       
tr.q      = 2;         % factor for radius change
tr.ip = @(x,y) d.dt*(x'*y); % inner product for CG iteration
tr.cgtol  = 1e-6;      % desired reduction of residual in CG
tr.cgits  = 50;        % maximum number of CG iterations

%% define target magnetization
% problem parameters
phase_shift = 0;             % alternating phase (pi) shifted excitation
sms         = 6;             % number of simultaneous slices (2,..,6)
slice_sep   = 0.025;         % multislice parameters
d.phi       = 90;            % flip angle in deg

% define center positions for all simultaneous slices 
if mod(sms,2) % (different for even/odd slice number)
    center_pos = slice_sep*((1:sms)-floor(sms/2+1))'; % odd
else
    center_pos = slice_sep*[fliplr(1:sms/2)-1/2 (-1)*((1:sms/2)-1/2)]'; % even
end

% create superposed slice profiles
inslice = sum(abs(repmat(d.xdis,size(center_pos,1),1)+...
                    repmat(center_pos,1,size(d.xdis,2)))<d.z);

% filter the target profile with a gaussian function
sigma       = 0.025;  % Gaussian filter paramters
filter_size = 75;
x = linspace(-filter_size,filter_size,length(inslice));
gaussFilter = exp(-x.^2/(2*sigma^2));
gaussFilter = gaussFilter/sum(gaussFilter); % normalize
d.inslice   = conv(double(inslice), gaussFilter, 'same');
d.outslice  = (1-d.inslice);

% create desired magnetization based on phi and inslice, outslice
d.Md = d.M0c*(repmat([0;sin(d.phi*pi/180);cos(d.phi*pi/180)],1,d.Nx).*repmat(d.inslice,3,1) + ...
              repmat([0;0;1],1,d.Nx).*(repmat(d.outslice,3,1))); 

% change desired magnetization by moving every second slice from y-plane 
% to x-plane which will add a phase alteration of pi 
if (phase_shift==1) && (sms>1)
    if mod(sms,2)
        mask = repmat((-1).^(0:sms-1)',1,size(d.xdis,2)); %odd
    else
        mask = repmat((-1).^(1:sms)',1,size(d.xdis,2)); %even
    end

    % create superposed slice profiles
    inslice = sum(mask.*(abs(repmat(d.xdis,size(center_pos,1),1)+...
                            repmat(center_pos,1,size(d.xdis,2)))<d.z) );
    
    % filter the target profile with a gaussian function
    d.inslice = conv(inslice, gaussFilter, 'same'); 

    % add a phase shift of pi/2 to even slice numbers (y->x) to preserve
    % a symmetric slice pattern; the additional phase has to be considered 
    % in the reconstruction 
    if mod(sms,2)
        d.Md(1,:) = 0; 
        d.Md(2,:) = d.inslice*sin(d.phi*pi/180); % odd slice number
    else
        d.Md(1,:) = d.inslice*sin(d.phi*pi/180); % even number          
        d.Md(2,:) = 0;  
    end
end

%% optimization
disp(['Computing minimizer for alpha = ' num2str(d.alpha)]);
u = tr_newton(d,tr,@objfun,@applyHess,d.u0);

%% plot results
plot_results(u,d);    
