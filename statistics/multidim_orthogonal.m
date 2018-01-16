function [ linco, offset,  Rsq , hist] = multidim_orthogonal( Y,X,nimp )
% [ linco, offset, res ] = multidim_orthognal( Y,X,nimp )
%   USES PCA to do a robus regression of Y from nimp columns in X
% 
%
% Y is the observable (vector N_experiments x 1)
% X are the parameters (matrix N_experiments x N_params)
% nimp is the number of predictors to be taken from N_params
%
% linco are the linear coefficients
% note that they are of the form ( Coeff1 Coeff2 Coeff3 ... )
% offset is the constant offset
% hist is the order of parameters used for regression
% Rsq is the Rsq value after N steps
% 
% Not suitable for anything, really
% Serge Dmitrieff
% www.biophysics.fr

S=size(Y);
if min(S)>1
    error('Now limited to 1 measured output')
else
    if S(1)==1
        Y=Y';
        S=S(2:-1:1);
    end
    ny=S(1);
end
S=size(X);
if S(2)==ny
    X=X';
    S=S(2:-1:1);
end
nx=S(2)+1;
[X,CX]=nono(X);
X=[Y X];


%% Preparing the data
chosable=ones(1,nx);
chosable(1)=0;
%factors=ones(1,nx);
%offsets=ones(1,nx);

for i=2:nx
	rangex=max(X(:,i))-min(X(:,i));
    if rangex==0
        chosable(i)=0;
		%else
    end
end


nvx=sum(chosable);
bli=1:nx;
%choices=bli(chosable);
        
if nimp>nvx
    nimp=nvx;
    disp('Warning : changing predictor number to variable number')
end       



left_ix=logical(chosable);
left=bli(left_ix);
hist=[1];
res=zeros(1,nimp);



%% Simple method
% Removing mean 

% We find the best to worse predictor in X
for ni=1:nimp
    scores=ones(1,nvx);
    
    
    for j=1:sum(left_ix)
        vec=[hist left(j)];
        M=X(:,vec);
        %coefs=M\YY;
		[~,~,latent]=princom(M);
        scores(j)=latent(end);
        %scores(j)=-sum(sc2);
    end
    
    [~,rk]=min(scores);
    ix=left(rk);
    % The next best predictor is #ix
    hist=[hist ix];
	res(ni)=scores(rk);
    % Updating variables we can still use
    left_ix(ix)=0;
    left_ix=logical(left_ix);
    left=bli(left_ix);
    nvx=sum(left_ix);
end
    
 M=X(:,vec);
[coefs]=princom(M);


[xprod]=ndim_xprod(coefs(:,1:ni));
linco=-(xprod(2:end)./CX(:))/xprod(1);

%coefs(:,1:ni)
	
%linco=coefs(2:end);
offset=coefs(1);
Rsq=1-res;
hist=hist(2:end)-1;

end

function [X,cof] = nono(X)
    s=size(X);
    b=ones(s(1),1);
    size(b*mean(X,1))
    X=X-b*mean(X,1);
    cof=var(X,1);
    X=X./(b*cof);
end
