function [m,J]=KProducts(z,K); 
% M.Terré, L. Féty, N. Paul
% Sept 2006
z=z(:).';

% Matrice des moments (en complexe)
U=ones(1,length(z));
for k=1:K-1
  U=[U;z.^k];
end
v=z.^K;

Rxx=U*U';
rxd=U*v';
a=Rxx\rxd;
ro=v*v';
J=real(ro-a'*rxd);
m=roots(fliplr([a' -1]));

