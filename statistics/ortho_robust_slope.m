function [ slope, off, score ] = ortho_robust_slope(PTS, verbose)
%function [ fit, score ] = ortho_robust_fit(PTS, verbose)
%
% give the slope of the line fitting the cloud of points PTS
% PTS is as vector:
%   x1 x2 x3 ... xn
%   y1 y2 y3 ... yn   (or in the other direction)
% Uses robust orthogonal fitting
% slope is the slope, off the offset
% 
% after F. Nédélec's ortho_robus_fit
% after Serge Dmitrieff's ortho_robust_coeff2
% S. Dmitrieff 10/01/2018
% www.biophysics.fr

if nargin < 2
    verbose = 0;
end

if size(PTS,1) == 2;
    PTS = PTS';
end

cen = mean(PTS, 1);

if numel(cen) ~= 2
    error('unexpected argument size');
end

    function res = residual(x)
        d = [cos(x(1)), sin(x(1))];
        a = cen + x(2) * [-d(2), d(1)];
        h = ( PTS(:,1) - a(1) )*d(2) - ( PTS(:,2) - a(2) )*d(1);
        % robust fitting:
        res = sum(abs(h));
        % for standard fitting, use this:
        % res = sum(h.^2);
    end

[coef, score] = fminsearch(@residual, [0,0]);

angle = coef(1);
offset = coef(2);

dir = [cos(angle), sin(angle)];
ori = cen + offset * [-dir(2), dir(1)];


slope=tan(angle);
off=ori(2)-ori(1)*slope;
%fit = [angle, ori(1), ori(2) ];

if verbose > 0
    
    figure('Position', [20 20 768 768]);
    plot(PTS(:,1), PTS(:,2), 'o');
    hold on;
    plot(cen(1), cen(2), 'ko', 'MarkerSize', 16);
    axis equal
    xlim([-1.0 4]);
    ylim([-3.5 1.5]);

    dif = max(PTS, [], 1) - min(PTS, [], 1);
    len = sqrt(sum(dif.^2)) / 2;
    
    % plot best line fit
    lin = vertcat(ori - len * dir, ori + len * dir);
    plot(lin(:,1), lin(:,2), '-', 'LineWidth', 2);
    
    % plot line fit going through the center:
    lin = vertcat(cen - len * dir, cen + len * dir);
    plot(lin(:,1), lin(:,2), '--', 'LineWidth', 2);
   
    % plot linear regression
    %poly = polyfit(PTS(:,1), PTS(:,2), 1);
    %val = [min(PTS(:,1)) max(PTS(:,1))];
    %plot(val, polyval(poly,val), ':', 'LineWidth', 1);

    % plot residual as a function of the angle
    val = 0:0.01:pi;
    res = zeros(size(val));
    for i = 1:length(val)
        res(i) = residual([val(i), 0]);
    end
    axes('Position', [0.06 0.825 0.2 0.15]);
    plot(val, res);
    xlim([0 pi]);
    xlabel('Angle (radian)');
    ylabel('Residual');

    % plot residual as a function of the offset
    val = -1:0.01:1;
    res = zeros(size(val));
    for i = 1:length(val)
        res(i) = residual([angle, val(i)]);
    end
    axes('Position', [0.06 0.625 0.2 0.15]);
    plot(val, res);
    xlabel('Perpendicular offset');
    ylabel('Residual');

    if verbose > 1

        % display orthogonal projections of all data points:
        a = ( PTS(:,1)-ori(1) )*dir(1) + ( PTS(:,2)-ori(2) )*dir(2);
        pX = ori(1) + a * dir(1);
        pY = ori(2) + a * dir(2);
        for n = 1 : length(PTS)
            line([pX(n), PTS(n,1)], [pY(n), PTS(n,2)]);
        end

    end
end


end



