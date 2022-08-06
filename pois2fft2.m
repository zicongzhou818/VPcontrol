function U=pois2fft2(F)
% h=1
% Solve 2-d Poisson equation
%
%      -div grad u = f on [0,1]x[0,1]
%
% with 0 boundary condition, i.e., u({0,1},y)=u(x,{0,1})=0.
%
% Method:
%
%   1) Set up the grid by partitioning along x- and y-directions
%      equidistant:
%
%         h=1/(n+1); xi=i*h; yj=j*h.
%
%   2) Let U(n-by-n) and F(n-by-n) with U(i,j)=u(xi,yj) to be determined
%      and F(i,j)=f(xi,yj). Then U satisfies
%
%            Tn*U+U*Tn=h^2*F,
%
%      where Tn(n-by-n)=tridiag(-1,2,-1).
%
%   3) Solve Tn*V+V*Tn=h^2*F by FFT (direct sine transformation)
%
%  Ref: http://www.ec-securehost.com/SIAM/ot56.html, section 6.7
%
%    RCL 11/15/2009
%
% Input
%
%     F   matrix (n-by-n)
%         F(i,j)=f(xi,yj)
%
% Output
%
%     U   matrix (n-by-n)
%         U(i,j) approximates u(xi,yj)

[m,n]=size(F);
if m~=n
   disp('pois2ftt2: F must be a square matrix');
   return
end
%h=1/(n+1); G=h^2*F;

U=dst(F.',n); U=U.'; U=dst(U,n);

theta=(pi/(2*(n+1)))*(1:n);
S=(4*(sin(theta')).^2)*ones(1,n)+ones(n,1)*(4*(sin(theta)).^2);
U=U./S;

U=idst(U.',n); U=U.'; U=idst(U,n);
