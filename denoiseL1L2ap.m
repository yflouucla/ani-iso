function u =denoiseL1L2ap(f, const, pm)
% |x|+|y|-alpha*sqrt(x^2+y^2) + 0.5*lambda*||Dx u-x+bx||^2 + 0.5*lambda*||Dy u-y+by||^2 + 0.5*mu*||u-f||^2
%


[rows,cols] = size(f);

mu = 20; lambda = 1; nIter = 200; u_orig = zeros(rows, cols);
u_orig = zeros(rows, cols); maxDCA = 10;
u0 = zeros(rows,cols); tol = 1e-4;

if isfield(pm,'lambda'); lambda = pm.lambda; end
if isfield(pm,'mu'); mu = pm.mu; end
if isfield(pm,'nIter'); nIter = pm.nIter; end
if isfield(pm,'maxDCA'); maxDCA = pm.maxDCA; end
if isfield(pm,'u_orig'); I = pm.u_orig; end
if isfield(pm,'u0'); u0 = pm.u0; end
if isfield(pm,'tol'); tol = pm.tol; end; % inner iteration tolerance




eps = 1e-16;
u = u0;
ux = Dx(u);
uy = Dy(u);
ugrad = sqrt(abs(ux).^2+abs(uy).^2);

oit = 1;
stop = 0;


% Build Kernel
uker = zeros(rows,cols);
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;uker(rows,1)=-1;uker(1,cols)=-1;
uker = mu+lambda*fft2(uker);

x = zeros(rows,cols);
y = zeros(rows,cols);


F(1) = sum(sum(abs(ux)+abs(uy)-const*ugrad)) + mu/2*norm(u-f,'fro')^2;

ff = f;
kkk = 1; tstart = tic;

while (oit <= maxDCA && stop == 0)
    uold2 = u;
    
    % Reserve memory for the auxillary variables
    bx = zeros(rows,cols);
    by = zeros(rows,cols);
    
    for  inner = 1:nIter
        uold = u;
        % update u
        rhs = mu*ff+lambda*Dxt(x-bx)+lambda*Dyt(y-by);
        u = real(ifft2(fft2(rhs)./uker));
        
        
        % update x and y
        dx = Dx(u);
        dy = Dy(u);
        % anisotropic TV
        x = shrink(dx+bx+const*ux/lambda./(ugrad+eps), 1/lambda);
        y = shrink(dy+by+const*uy/lambda./(ugrad+eps), 1/lambda);
        
        % update bregman parameters
        bx = bx+dx-x;
        by = by+dy-y;
        
        if (norm(u-uold, 'fro')/norm(uold,'fro')<tol)
            break;
        end
        
    end
    
    
    output.relerr(oit)=norm(u-uold2,'fro');
    
    ux = Dx(u);
    uy = Dy(u);
    
    ugrad = sqrt(abs(ux).^2+abs(uy).^2);
    F(oit+1) = sum(sum(abs(ux)+abs(uy)-const*ugrad)) + mu/2*norm(u-f,'fro')^2;
    
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

function [xs,ys] = shrink2(x,y,lambda)

s = sqrt(x.*conj(x)+y.*conj(y));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;

return;

function z = shrink(x,r)
z = sign(x).*max(abs(x)-r,0);
return;

