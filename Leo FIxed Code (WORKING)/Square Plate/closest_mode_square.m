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

rows = [];
for f = f_meas(:).'
    k2drive = 2*pi*f/cf;                        % k^2 at the drive frequency
    Nmax    = ceil(L/pi*sqrt(2*k2drive)) + 2;   % search up to twice the drive k^2
    [n, m]  = meshgrid(0:Nmax);
    driven  = abs(cos(n*pi*x0/L).*cos(m*pi*y0/L)) > 1e-9;
    keep    = (n >= m) & (n + m > 0) & driven;
    n = n(keep);  m = m(keep);

    fpred = cf*(pi/L)^2*(n.^2 + m.^2)/(2*pi);
    err   = 100*(fpred - f)/f;
    [~, idx] = sort(abs(err));
    idx = idx(1:min(nShow, numel(idx)));

    rows = [rows; repmat(f, numel(idx), 1), n(idx), m(idx), fpred(idx), err(idx)]; %#ok<AGROW>
end

T = array2table(rows, 'VariableNames', {'f_meas_Hz', 'n', 'm', 'f_pred_Hz', 'error_pct'});
end
