function [mv, sv, mo, so, slopes,offsets] = bootstrap_robust(PTS, n_max, verbose)
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
    verbose = 0;
end


if size(PTS,1) == 2;
    % This fails if PTS is 2x2.
    % If you do a bootstrap on 2x2 array, you deserve failure anyway
    %         (much more than if your data has the wrong orientation)
    PTS = PTS';
end

i_max = size(PTS, 1);
slopes = zeros(1, n_max);
offsets = zeros(1, n_max);
for n = 1 : n_max
    i = randi(i_max, i_max, 1);
    % Slower but more precise :
    [slopes(n),offsets(n)] = ortho_robust_slope(PTS(i, :));
    % Faster but worse for small number of pts :
    %[slopes(n),offsets(n)] = ortho_robust_coeff2(PTS(i, :)');
    % For comparison : linear fit
    %p=polyfit(PTS(i,1),PTS(i,2),1);
    %slopes(n)=p(1);
    %offsets(n)=p(2);
end

mv = mean(slopes);
sv = std(slopes);

mo = mean(offsets);
so = std(offsets);

if verbose
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
    
    figure
    scatter(PTS(:,1),PTS(:,2))
    hold all
    plot([min(PTS(:,1)) max(PTS(:,1))],[min(PTS(:,1)) max(PTS(:,1))]*mv+mo,'b')
end

end