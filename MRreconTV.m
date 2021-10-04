function u = MRreconTV(R,f, pm)
%function [u, err, cpu] = MRreconTV(R,f, mu, lambda,  nIter, ifcon, u_orig)

% |x|+|y| + 0.5*lambda*||Dx u-x+bx||^2 + 0.5*lambda*||Dy u-y+by||^2
% s.t. RFu = f
%
% using the Split Bregman method ftp://ftp.math.ucla.edu/pub/camreport/cam08-29.pdf
%


[rows,cols] = size(f);

mu = 20; lambda = 1; nIter = 200; u_orig = zeros(rows, cols);
u_orig = zeros(rows, cols);
u0 = zeros(rows,cols); tol = 1e-4;

if isfield(pm,'lambda'); lambda = pm.lambda; end
if isfield(pm,'mu'); mu = pm.mu; end
if isfield(pm,'maxit'); maxit = pm.maxit; end
if isfield(pm,'u_orig'); u_orig = pm.u_orig; end
if isfield(pm,'u0'); u0 = pm.u0; end
if isfield(pm,'tol'); tol = pm.tol; end; % inner iteration tolerance

[rows,cols] = size(f);

% Reserve memory for the auxillary variables
f0 = f;
u = u0;
x = zeros(rows,cols);
y = zeros(rows,cols);
bx = zeros(rows,cols);
by = zeros(rows,cols);

% Build Kernels
scale = sqrt(rows*cols);
murf = ifft2(mu*(conj(R).*f))*scale;

uker = zeros(rows,cols);
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
uker = mu*(conj(R).*R)+lambda*fft2(uker);


%  Do the reconstruction
for outer = 1:50
    for inner = 1:maxit/50       
        % update u
        rhs = murf+lambda*Dxt(x-bx)+lambda*Dyt(y-by);
        u = ifft2(fft2(rhs)./uker);
                
        % update x and y
        dx = Dx(u);
        dy  =Dy(u);
        
        % anisotropic TV
        x = shrink(dx+bx, 1/lambda);
        y = shrink(dy+by, 1/lambda);
        
        % update bregman parameters
        bx = bx+dx-x;
        by = by+dy-y;
    end
    
    f = f+f0-R.*fft2(u)/scale;
    murf = ifft2(mu*R.*f)*scale;
    
end


return;


function d = Dx(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);
d(:,1) = u(:,1)-u(:,cols);
return

function d = Dxt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);
d(:,cols) = u(:,cols)-u(:,1);
return

function d = Dy(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);
d(1,:) = u(1,:)-u(rows,:);
return

function d = Dyt(u)
[rows,cols] = size(u);
d = zeros(rows,cols);
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);
d(rows,:) = u(rows,:)-u(1,:);
return


function z = shrink(x,r)
z = sign(x).*max(abs(x)-r,0);
return;

