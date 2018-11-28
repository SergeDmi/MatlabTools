function [ y_fit , y_std , x_new ] = rolling_average( Y,X,dist,x_fit)
%ROLLING_AVERAGE Performs a rolling average
%   Gives y_fit corresponding to bins x_fit from data Y,X 
%   with a haracteristic distance dist
%   for now only gaussian weight
%
% INPUT
%   X is a (nf * 1) vector of coordinates that will be rolled on
%   Y is an (nf * ny ) array of observables to be be averaged
%   dist is the rolling distance
%   xfit is the fitted coordinates ; takes X if empty.
%
% OUTPUT
%   y_fit is the averaged Y
%   y_std is the standard deviation
%   x_new is the actual vector of x from the moving average
%
% @TODO : add function handle as input parameter
%
% Serge Dmitrieff, Institut Jacques Monod, 2018
% www.biophysics.fr

if nargin<4
	x_fit=X;
end

% First we put data in the correct format
sX=length(X);
nf=length(x_fit);
x_fit=reshape(x_fit,1,nf);
sY=size(Y);
X=reshape(X,sX,1);
if sY(2)==sX
	Y=Y';
	sY=size(Y);
end

nY=sY(2);
y_fit=zeros(nf,nY);
y_fit=zeros(nf, 1);
y_std=zeros(nf,nY);

ro=ones(1,nf);
XX=X*ro;
col=ones(sX,1);
weights=exp(-((XX-col*x_fit).^2)./(2*dist*dist));

x_new(:)=sum( (XX.*weights),1 )./sum( weights,1);

for i=1:nY
	% Making matrixes rather than loops
	YY=Y(:,i)*ro;
	% 3,2,1, let's jam
	y_fit(:,i)=sum( (YY.*weights),1 )./sum( weights,1);
	y_std(:,i)=sqrt(sum( ((YY-col*(y_fit(:,i)')).^2).*weights,1)./sum(weights,1)); 
end
% That was fast !

end

