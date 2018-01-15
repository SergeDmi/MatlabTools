function [xprod]=ndim_xprod(M)
%    [xprod]=ndim_xprod(M)
% A code to do an n-dimensional cross-product of n-1 vectors
%
% Serge Dmitrieff, IJM 2018
% www.biophysics.fr

sm=size(M);
if sm(2)==sm(1)-1
    M=M';
    N=sm(1);
elseif sm(1)==sm(2)-1
    %We good
    N=sm(2);
else
    error('Input matrix size should be N x (N-1) ');
end

d=1;
% we try to add a vector 
xprod=get_xprod(M,N,d);
while sum(isnan(xprod)) && d<N
    d=d+1;
    xprod=get_xprod(M,N,d);
end

if sum(isnan(xprod))
    error('Invalid input matrix : cannot find adjoint matrix')
end

end

function [xp]=get_xprod(M,N,d)
    MM=[M ; zeros(1,N)];
    MM(N,d)=1;
    % Computing adjoint matrix
    AA=det(MM)*inv(MM);
    xp=AA(:,end);
end