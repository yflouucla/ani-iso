function u = MRreconL1L2ap(R,f, pm)
% |x|+|y| - alpha sqrt(x^2+y^2) + 0.5*lambda*||Dx u-x+bx||^2 + 0.5*lambda*||Dy u-y+by||^2
% s.t. RFu = f
%
% DCA and Split Bregman
%

[rows,cols] = size(f);

mu = 20; lambda = 5; alpha = 0.5;
maxit = 1000;
u_orig = zeros(rows, cols);
u0 = zeros(rows,cols); tol = 1e-4;

if isfield(pm,'mu'); mu = pm.mu; end
if isfield(pm,'lambda'); lambda = pm.lambda; end
if isfield(pm,'alpha'); alpha = pm.alpha; end
if isfield(pm,'maxit'); maxit = pm.maxit; end
if isfield(pm,'u_orig'); u_orig = pm.u_orig; end
if isfield(pm,'u0'); u0 = pm.u0; end
if isfield(pm,'tol'); tol = pm.tol; end; % inner iteration tolerance


u = u0;
ux = Dx(u);
uy  =Dy(u);
ugrad = sqrt(abs(ux).^2+abs(uy).^2);
eps = 1e-8;

%DCA iteration counts
oit = 1;
stop = 0;
kkk = 1;
tstart = tic;


% Build Kernels
scale = sqrt(rows*cols);
murf = ifft2(mu*(conj(R).*f))*scale;

uker = zeros(rows,cols);
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
uker = mu*(conj(R).*R)+lambda*fft2(uker);

x = zeros(rows,cols);
y = zeros(rows,cols);

ff = f;

F(1) = sum(sum(abs(ux)+abs(uy)-alpha*ugrad)) + mu/2*norm(R.*fft2(u)/scale-f)^2;



while (oit <= 10 && stop == 0)
    
    uold = u;
    
    % Reserve memory for the auxillary variables
    f0 = ff; f = f0; murf = ifft2(mu*R.*f0)*scale;
    
    bx = zeros(rows,cols);
    by = zeros(rows,cols);
    
    
    
    %  Do the reconstruction
    for outer = 1:50
        for inner = 1:maxit/50
            % update u
            rhs = murf+lambda*Dxt(x-bx)+lambda*Dyt(y-by);
            u = (ifft2(fft2(rhs)./uker));
            
            % update x and y
            dx = Dx(u);
            dy = Dy(u);
            
            % anisotropic TV
            x = shrink(dx+bx+alpha*ux/lambda./(ugrad+eps), 1/lambda);
            y = shrink(dy+by+alpha*uy/lambda./(ugrad+eps), 1/lambda);
            
            % update bregman parameters
            bx = bx+dx-x;
            by = by+dy-y;
            
        end
        
        f = f+f0-R.*fft2(u)/scale;
        murf = ifft2(mu*R.*f)*scale;
        
    end
    
    u = abs(u);
    ux = Dx(u);
    uy = Dy(u);
    ugrad = sqrt(abs(ux).^2+abs(uy).^2);
    
    err2(oit) = norm(uold-u,'fro');
    
    F(oit+1) = sum(sum(abs(ux)+abs(uy)-alpha*ugrad))+ mu/2*norm(R.*fft2(u)/scale-f)^2;
    
    
    oit = oit + 1;
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

