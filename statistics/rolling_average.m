function [ yfit ] = rolling_average( Y,X,xfit,dist )
%ROLLING_AVERAGE Performs a rolling average
%   Gives yfit corresponding to bins xfit from data Y,X 
%   with a haracteristic distance dist
%   for now only gaussian weight
%
% Serge Dmitrieff, Institut Jacques Monod, 2018
% www.biophysics.fr

% First we put data in the correct format
sX=length(X);
nf=length(xfit);
xfit=reshape(xfit,1,nf);
X=reshape(X,sX,1);
Y=reshape(Y,sX,1);

% Making matrixes rather than loops
ro=ones(1,nf);
XX=X*ro;
YY=Y*ro;

% 3,2,1, let's jam
weights=exp(-((XX-ones(sX,1)*xfit).^2)./(2*dist*dist));
yfit=sum( (YY.*weights),1 )./sum( weights,1);

% That was fast !

end

