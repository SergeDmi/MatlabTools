function [MV, SV, MO, SO, Rsq,hist] = bootstrap_multi(Y,X, nimp, n_max, verbose)
%function [MV, SV, MO, SO, Rsq,hist] = bootstrap_multi(PTS, n_max)
%
% gives the mean (MV) and std dev (SV) of linear coeffs obtained from PTS
% gives the mean (MO) and std dev (SO) of y-offsets obtained from PTS
% Slopes & offsets are calculated using least-square fitting (sorry no robust)
% Y is the observable (vector N_experiments x 1)
% X are the parameters (matrix N_experiments x N_params)      
% n_max is the number of bootstraped samples. Higher is better.
% nimp is the number of predictors to be taken from N_params
%
% Based on F. Nédélec's Bootstrap implementation for a orthogonal fitting
% Serge Dmitrieff - 10/01/2018
% www.biophysics.fr

if nargin < 5
    verbose = 0;
end

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
nx=S(2);


slopes = zeros(nx, n_max);
offsets = zeros(1, n_max);
counts=zeros(nx,1);
r2s=zeros(nx,1);
bli=1:nx;

%% The actually bootstrap : reshuffling data
for n = 1 : n_max
    i = randi(ny, ny, 1);
    [sl,offsets(n), Rsq , hist] = multidim_regression( Y(i),X(i,:),nimp );
    slopes(hist,n)=sl;
    counts(hist)=counts(hist)+1;
    r2s(hist)=r2s(hist)+diff([0 Rsq])';
end
r2s=r2s./counts;

%% Getting means and variance
MV=zeros(1,nx);
SV=zeros(1,nx);
for x=1:nx
    slope=slopes(x,:);
    sl=slope(slope~=0);
    % We count in the average only if is actually a computed value, i.e. ~=0
    MV(x)=mean(sl);
    SV(x)=std(sl);
end
MO = mean(offsets);
SO = std(offsets);

%% Plotting if verbose
if verbose
    figure
    bar(bli,counts);
    xlabel('component')
    ylabel('Times used in regression')
    figure
    bar(bli,r2s)
    ylabel('Average R^2')
    xlabel('component')
    for x=1:nx
        figure
        title(['component #' num2str(x)])
        h=histogram(slopes(x,:), 32);
        hm=max(h.Values);
        xlabel('Coefficient');
        ylabel('Histogram count');
        hold all
        % mean slope +/- stdev
        plot([MV(x) MV(x)],[0 hm],'b','LineWidth',3)
        plot([MV(x)-SV(x) MV(x)-SV(x)],[0 hm],'r','LineWidth',2)
        plot([MV(x)+SV(x) MV(x)+SV(x)],[0 hm],'r','LineWidth',2)
    end
end

end