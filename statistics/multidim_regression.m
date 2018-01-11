function [ linco, offset,  Rsq , hist] = multidim_regression( Y,X,nimp )
% [ linco, offset, res ] = multidim_regression( Y,X,nimp )
%   A function to do a linear regression of Y from nimp columns in X
% 
%
% Y is the observable (vector N_experiments x 1)
% X are the parameters (matrix N_experiments x N_params)
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
X=[ones(ny,1) X];


%% Preparing the data
chosable=ones(1,nx);

%factors=ones(1,nx);
%offsets=ones(1,nx);

for i=1:nx
	rangex=max(X(:,i))-min(X(:,i));
    if rangex==0
        chosable(i)=0;
    else
        if 0
            % this is if you want to renormalize your data
            % warning, you should also renormalize slopes afterwards
            factors(i)=rangex;
            offsets(1,i)=mean(X(:,i));
            X(:,i)=X(:,i)-offsets(1,i);
            X(:,i)=X(:,i)/rangex;
        end
    end
end

nvx=sum(chosable);
bli=1:nx;
%choices=bli(chosable);
        
if nimp>nvx
    nimp=nvx;
    disp('Warning : change predictor number to variable number')
end       



left_ix=logical(chosable);
left=bli(left_ix);
hist=zeros(1,nimp);
res=zeros(1,nimp);



%% Simple method
% Removing mean 
M=X(:,1);
cte=M\Y;
YY=Y-M*cte;
err0=sum(YY.^2);
% We find the best to worse predictor in X
for ni=1:nimp
    scores=ones(1,nvx);
    coeffs=ones(2,nvx);
    
    for j=1:sum(left_ix)
        vec=[1 left(j)];
        M=X(:,vec);
        coefs=M\YY;
        dy=YY-M*coefs;
        scores(j)=sum(dy.^2);
        coeffs(:,j)=coefs;
    end
    
    [~,rk]=min(scores);
    ix=left(rk);
    % The next best predictor is #ix
    hist(ni)=ix;
    % We recompute the fitting error 
    % This should result in better fitting but is slower
    vec=[1 hist(1:ni)];
    M=X(:,vec);
    coefs=M\Y;
    YY=Y-M*coefs;
    res(ni)=sum(YY.^2);
    % Updating variables we can still use
    left_ix(ix)=0;
    left_ix=logical(left_ix);
    left=bli(left_ix);
    nvx=sum(left_ix);
end
    
linco=coefs(2:end);
offset=coefs(1);
Rsq=1-res/err0;
hist=hist-1;

end

