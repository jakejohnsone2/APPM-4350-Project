function T = closest_mode_circle(f_meas, cf, R, nShow)
%CLOSEST_MODE_CIRCLE  Nearest circular-plate resonances that a center force can drive.
%   T = closest_mode_circle(f_meas, cf, R)          nearest ring mode for each f
%   T = closest_mode_circle(f_meas, cf, R, nShow)   nShow candidates per f
%
%   f_meas : measured frequency or vector of frequencies [Hz]
%   cf     : sqrt(D/(rho*h)), so f = cf*k^2/(2*pi)
%   R      : plate radius [m]
%
% Model: f_ns = cf*(j'_ns/R)^2/(2*pi), where j'_ns is the s-th positive zero of J_n'.
% A center force only drives modes that are nonzero at r = 0. Since J_n(0) = 0 for
% every n >= 1, only the n = 0 (ring) modes are listed. The rigid mode (k = 0) is
% excluded. Needs circ_roots.m.

if nargin < 4 || isempty(nShow), nShow = 1; end

% Pre allocating the matrix
rows = [];

% For loop that goes through each of our lab frequencies
for f = f_meas(:).'
    % Finding the value of the derivative of the Bessel function at our
    % frequency
    xneed = R*sqrt(2*pi*f/cf);

    % Setting the max for the roots were going to search
    smax  = ceil(2*xneed/pi) + 3;            

    % Calculating the zeros of our Bessel Function
    jp    = circ_roots(0, smax);     

    % Setting the index of each root
    s     = (1:smax).';

    % Predicting the frequency of our surviving modes
    fpred = cf*(jp(1, :).'/R).^2/(2*pi);

    % Finding the error between the predicted frequency and lab frequency
    err   = 100*(fpred - f)/f;


    % Sorting by the smallest error
    [~, idx] = sort(abs(err));
    idx = idx(1:min(nShow, numel(idx)));
    
    % Appending the matrix such that we are only left with one candidate
    % per lab frequency
    rows = [rows; repmat(f, numel(idx), 1), zeros(numel(idx), 1), s(idx), fpred(idx), err(idx)]; %#ok<AGROW>
end

% Outputing the final table
T = array2table(rows, 'VariableNames', {'f_meas_Hz', 'n', 's', 'f_pred_Hz', 'error_pct'});
end
