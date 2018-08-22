function [ lincoeffs, offset,  Rsq , hist] = multidim_orthogonal( Y,X,nimp )
% [ lincoeffs, offset, res ] = multidim_orthognal( Y,X,nimp )
%   USES PCA to do a robus regression of Y from nimp columns in X
% 
%
% Y is the observable (vector N_experiments x 1)
% X are the parameters (matrix N_experiments x N_params)
% nimp is the number of predictors to be taken from N_params
%
% lincoeffs are the linear coefficients
% note that they are of the form ( Coeff1 Coeff2 Coeff3 ... )
% offset is the constant offset
% hist is the order of parameters used for regression
% Rsq is the Rsq value after N steps
% 
% Not suitable for anything, really
% Serge Dmitrieff 2018
% www.biophysics.fr

%% Checking input.
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
XX=[Y X];
[XX,CX]=nono(XX);

varY=var(Y);

%% Preparing the data
chosable=ones(1,nx);
chosable(1)=0;
chosable(2:end)=CX(2:end)>0;

nvx=sum(chosable);
bli=1:nx;
        
if nimp>nvx
    nimp=nvx;
    disp('Warning : changing predictor number to variable number')
end       

left_ix=logical(chosable);
left=bli(left_ix);
hist=[1];
res=zeros(1,nimp);

%% Simple method

% We find important variables in order
for ni=1:nimp
    scores=ones(1,nvx);
    
    % We find the best to worse predictor in X
    for j=1:sum(left_ix)
        vec=[hist left(j)];
        M=XX(:,vec);
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
    

 M=XX(:,hist);
 % Back to computing CPA
[coefs,PP,latent]=princom(M);
% Computing normal vector to main plane and getting linear coeffcients
[xprod]=ndim_xprod(coefs(:,1:ni));
linco=-xprod(2:end)/xprod(1);
% Back to real coefficients by inverting normalization
hist=(hist(2:end)-1);
linco=CX(1)*linco'./CX(hist+1);
lincoeffs=zeros(1,nx-1);
lincoeffs(hist)=linco;
%normdv=xprod./sqrt(sum(xprod(:).^2));
%normdv=xprod;
%err=latent(end)*sqrt(sum(xprod.*(CX(1:(ni+1)).^2)'));
%normdv=normdv.*CX(1:(ni+1))';
%allvecs=PP(:,end)*normdv';
%pp=PP;
%pp(:,1)=0;
%ppp=pp*inv(coefs);

%allvecs(:,1)=ppp(:,2)*CX(2);
%allvecs(:,2)=ppp(:,1)*CX(1);

predY=sum((ones(ny,1)*(lincoeffs)).*X,2);
predY=predY-mean(predY)+mean(Y);
errs=Y-predY;
offset=mean(errs);


%Rsq(nimp)=1-(devY-std(errs))/devY;
% Computing the Rsquared
Rsq=(varY-var(errs))/varY;

%scatter(XX(:,2)*CX(2),XX(:,1)*CX(1));
%rest=mean(Y)-sum(mean(X,1).*linco);
%offset=coefs(1);


end

function [X,cof] = nono(X)
    % Normalizing input data :
    % centering on zero and dividing by deviation
    s=size(X);
    b=ones(s(1),1);
    X=X-b*mean(X,1);
    cof=std(X,1);
    X=X./(b*cof);
end
