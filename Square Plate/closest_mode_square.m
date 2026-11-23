function T = closest_mode(f, f11, nShow, doPlot)
%CLOSEST_MODE  Nearest model resonances (n,m) to a measured frequency.
%   T = closest_mode(f)                 f in Hz, calibrated with (1,1) = 175 Hz
%   T = closest_mode(f, f11)            calibrate with your measured (1,1) [Hz]
%   T = closest_mode(f, f11, nShow)     number of candidates to list (default 5)
%   T = closest_mode(f, f11, nShow, 0)  skip the plots
%
% Model: f = C*(n^2 + m^2), with C = f11/2 (from the square-plate Helmholtz
% model, calibrated on the (1,1) mode). Candidates are listed by smallest
% percent error. Modes with the same n^2 + m^2 (e.g. (5,0) and (4,3)) are
% degenerate in the model, so they show up with identical predicted frequencies.
%
% center_node = true when n or m is odd, i.e. W(center) = 0. These are the modes
% a plate pinned at the center by a screw can resonate in.
%
% Plots: the symmetric (+) and antisymmetric (-) combination for the top
% candidates, W = cos(n*pi*x/L)cos(m*pi*y/L) +/- cos(m*pi*x/L)cos(n*pi*y/L).
%
% Example:  closest_mode(1608)      closest_mode(3201, 175, 6)

if nargin < 2 || isempty(f11),   f11   = 175; end
if nargin < 3 || isempty(nShow), nShow = 5;   end
if nargin < 4 || isempty(doPlot), doPlot = true; end

L  = 0.1334;                         % plate side [m]
rc = 0.003;                          % screw head radius [m] (masked in plots)
C  = f11/2;                          % [Hz] per unit of n^2 + m^2

% Candidate modes (n >= m, excluding the rigid-body mode)
Nmax = ceil(sqrt(2*f/C)) + 2;
[n, m] = meshgrid(0:Nmax, 0:Nmax);
keep = (n >= m) & (n + m > 0);
n = n(keep);  m = m(keep);

fpred = C*(n.^2 + m.^2);
err   = 100*(fpred - f)/f;
[~, idx] = sort(abs(err));
idx = idx(1:min(nShow, numel(idx)));

n_col = n(idx);  m_col = m(idx);
T = table(n_col, m_col, fpred(idx), err(idx), mod(n_col,2)==1 | mod(m_col,2)==1, ...
          'VariableNames', {'n', 'm', 'f_pred_Hz', 'error_pct', 'center_node'});
fprintf('Measured f = %.0f Hz   (C = %.2f Hz, from (1,1) = %.0f Hz)\n', f, C, f11);
disp(T);

if doPlot
    nP = min(3, height(T));
    x = linspace(0, L, 400);
    [X, Y] = meshgrid(x);
    figure('Position', [100 100 800 300*nP]);
    for i = 1:nP
        for j = 1:2
            s = 3 - 2*j;                           % j = 1 -> +1, j = 2 -> -1
            subplot(nP, 2, 2*(i-1) + j);
            ni = T.n(i);  mi = T.m(i);
            if s < 0 && ni == mi
                text(0.5, 0.5, 'n = m: no antisymmetric partner', ...
                     'HorizontalAlignment', 'center');
                axis off;
                continue
            end
            W = cos(ni*pi*X/L).*cos(mi*pi*Y/L) + s*cos(mi*pi*X/L).*cos(ni*pi*Y/L);
            Wabs = abs(W)/max(abs(W(:)));
            Wabs(hypot(X - L/2, Y - L/2) < rc) = NaN;   % hide the screw region
            imagesc(x, x, Wabs, 'AlphaData', ~isnan(Wabs)); axis image xy;
            colormap(flipud(gray)); caxis([0 0.15]);
            if s > 0, sg = '+'; else, sg = '-'; end
            title(sprintf('(%d,%d)%s, f_{pred} = %.0f Hz (%+.1f%%)', ...
                  ni, mi, sg, T.f_pred_Hz(i), T.error_pct(i)));
        end
    end
end
end
