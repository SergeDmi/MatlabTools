function [mv, sv, mo, so, mr,sr,slopes,offsets,rsqs] = bootstrap_robust(Y,X, n_max, verbose)
%function [mv, sv,slo] = bootstrap_robust(PTS, n_max)
%
% gives the mean (mv) and std dev (sv) of slopes obtained from PTS
% gives the mean (mo) and std dev (so) of y-offset obtained from PTS
% Slopes & offsets are calculated using orthogonal fitting
% PTS is as vector 
%   x1 x2 x3 ... xn
%   y1 y2 y3 ... yn        (or the orgthogonal)
% n_max is the number of bootstraped samples. Higher is better.
%
% Based on F. Nédélec's Bootstrap implementation
% Serge Dmitrieff - 10/01/2018
% www.biophysics.fr

if nargin < 3
	n_max=max(size(Y))*10;
end

if nargin < 4
    verbose = 0;
end

if size(X,1) == 1
    X = X';
end

if size(Y,1) == 1
    Y = Y';
end


i_max = size(Y, 1);
slopes = zeros(1, n_max);
offsets = zeros(1, n_max);
rsqs = zeros(1, n_max);
for n = 1 : n_max
    i = randi(i_max, i_max, 1);
    [slopes(n),offsets(n),~ ,rsqs(n)] = ortho_robust_slope(Y(i),X(i));
end

mv = mean(slopes);
sv = std(slopes);

mo = mean(offsets);
so = std(offsets);

mr=mean(rsqs);
sr=std(rsqs);

if verbose>0
	figure
	scatter(X,Y)
	hold all
	plot([min(X) max(X)],[min(X) max(X)]*mv+mo,'r')
	
    if verbose>1
        figure
		% Histogram of slope values
		h=histogram(slopes, 32);
		hm=max(h.Values);
		xlabel('Slope');
		ylabel('Histogram count');
		hold all
		% mean slope +/- stdev
		plot([mv mv],[0 hm],'b','LineWidth',3)
		plot([mv-sv mv-sv],[0 hm],'r','LineWidth',2)
		plot([mv+sv mv+sv],[0 hm],'r','LineWidth',2)
    end
  
end

end
