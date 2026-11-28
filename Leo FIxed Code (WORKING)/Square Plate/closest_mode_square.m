function T = closest_mode_square(f_meas, L, cf, x0, y0, nShow)
%CLOSEST_MODE_SQUARE  Nearest square-plate resonances that a force at (x0,y0) can drive.
%   T = closest_mode_square(f_meas, L, cf, x0, y0)          nearest mode for each f
%   T = closest_mode_square(f_meas, L, cf, x0, y0, nShow)   nShow candidates per f
%
%   f_meas : measured frequency or vector of frequencies [Hz]
%   L      : side length [m]
%   cf     : sqrt(D/(rho*h)), so f = cf*k^2/(2*pi)
%   x0, y0 : drive point [m] (L/2, L/2 for the center screw)
%
% Model: f_nm = cf*(pi/L)^2*(n^2 + m^2)/(2*pi).
% A mode is driven only if phi_nm(x0,y0) = cos(n*pi*x0/L)*cos(m*pi*y0/L) ~= 0,
% because the projection coefficient is proportional to that value (Report Eq. 9).
% At the center this keeps only n and m both even. Degenerate pairs (n,m) and (m,n)
% have the same frequency and are listed once, with n >= m. The rigid mode (0,0)
% is excluded.

if nargin < 6 || isempty(nShow), nShow = 1; end

% Pre allocating the matrix
rows = [];

% For loop that goes through each of our lab frequencies
for f = f_meas(:).'
    % Calculating k squared at the respective frequency
    k2drive = 2*pi*f/cf;                      

    % Setting the max n for how far we are going to search
    Nmax    = ceil(L/pi*sqrt(2*k2drive)) + 2; 

    % Creaing a meshgrid for the ns and ms
    [n, m]  = meshgrid(0:Nmax);

    % Finding the modes value at the origin for the respective frequency.
    % Filters out the odd numbers as for even numbers the value will be 0.
    driven  = abs(cos(n*pi*x0/L).*cos(m*pi*y0/L)) > 1e-9;
    
    % Filtering out such that we are left with only one of each pair, no
    % (0,0) mode, and the modes that cant be excited (odd modes)
    keep    = (n >= m) & (n + m > 0) & driven;

    % Reducing so we are only left with the filtered values
    n = n(keep);  m = m(keep);
    
    % Predicting the frequency of our surviving modes
    fpred = cf*(pi/L)^2*(n.^2 + m.^2)/(2*pi);

    % Finding the error between the predicted frequency and lab frequency
    err   = 100*(fpred - f)/f;

    % Sorting by the smallest error
    [~, idx] = sort(abs(err));
    idx = idx(1:min(nShow, numel(idx)));
    
    % Appending the matrix such that we are only left with one candidate
    % per lab frequency
    rows = [rows; repmat(f, numel(idx), 1), n(idx), m(idx), fpred(idx), err(idx)]; %#ok<AGROW>
end

% Outputing the final table
T = array2table(rows, 'VariableNames', {'f_meas_Hz', 'n', 'm', 'f_pred_Hz', 'error_pct'});
end
