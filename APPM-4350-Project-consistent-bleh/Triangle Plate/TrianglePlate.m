%% TrianglePlate.m -- Chladni figures on a free-edge equilateral triangle, center drive
% Model (Report Sec. 2.3, 2.6):
%   lap(W) + (k + i*gamma)^2 W = -delta(x-x0) delta(y-y0),   dW/dn = 0 on the edges
% with Lame's modes: each mode is a sum of 6 plane waves (tri_mode.m),
%   k_mn^2 = (4*pi/(3L))^2 (m^2 + m*n + n^2).
%
% The center screw is modeled as a point force at the centroid (averaged over the
% screw head). A mode is only driven if phi_mn(x0,y0) ~= 0, which at the centroid
% means m - n is a multiple of 3. Every pattern below is computed at a MEASURED lab
% frequency and labeled with the nearest mode the center force can drive.
%
% Light regions in the plot = |W| ~ 0 = where the sand collects.
% Needs tri_mode.m and closest_mode_triangle.m in the same folder.
clear; clc; close all

%% 1. Plate properties and frequency -> wave number (Sec. 2.2)
E   = 69e9;  nu = 0.33;  rho = 2700;  h = 0.52e-3;   % aluminum (thinner than the others)
H   = 0.1268;                                        % triangle height [m]
L   = 2*H/sqrt(3);                                   % side length [m]
D   = E*h^3/(12*(1 - nu^2));                         % bending stiffness
cf  = sqrt(D/(rho*h));                               % omega = cf*k^2

%% 2. Numerical settings (Sec. 2.7)
P.Ng  = 1000;               % grid points per side
P.gam = 0.01*pi/L;          % damping gamma [1/m]; keeps resonance finite
P.x0  = L/2;  P.y0 = H/3;   % center drive (centroid)
P.rc  = 0.0026;             % screw head radius [m]

%% 3. Measured frequencies from the lab [Hz]
f_lab = [230 600 1950 2250 3390 4835];

%% 4. Pattern at each measured frequency
for f = f_lab
    [Wabs, x, y, inside] = triangle_response(f, L, cf, P);
    M = closest_mode_triangle(f, L, cf, P.x0, P.y0, 1);   % nearest driven mode
    lbl = sprintf('f = %.0f Hz (nearest driven mode (%d,%d), %.0f Hz)', ...
                  f, M.m, M.n, M.f_pred_Hz);
    xt = [0, L, L/2, 0];  yt = [0, 0, H, 0];             % triangle outline

    % Figure 1: full mirrored pattern
    figure;
    imagesc(x, y, Wabs, 'AlphaData', ~isnan(Wabs)); axis image xy;
    colormap(flipud(gray)); caxis([0 0.15]);
    hold on; plot(xt, yt, 'r-', 'LineWidth', 1.5); hold off;
    title(['Triangle (mirrored), ' lbl]);
    xlabel('x [m]'); ylabel('y [m]');
    saveas(gcf, sprintf('triangle_full_pattern_f%dHz.png', round(f)));

    % Figure 2: just the plate
    Wtri = Wabs;  Wtri(~inside) = NaN;
    figure;
    imagesc(x, y, Wtri, 'AlphaData', ~isnan(Wtri)); axis image xy;
    colormap(flipud(gray)); caxis([0 0.15]);
    xlim([0 L]); ylim([0 H]);
    title(['Triangle, ' lbl]);
    xlabel('x [m]'); ylabel('y [m]');
    saveas(gcf, sprintf('triangle_pattern_f%dHz.png', round(f)));
end

%% 5. Table for the report (Table 1)
T = closest_mode_triangle(f_lab, L, cf, P.x0, P.y0, 1);
disp(T);

%% ------------------------------------------------------------------------
function [Wabs, x, y, inside] = triangle_response(f, L, cf, P)
%TRIANGLE_RESPONSE  Normalized |W| for a centroid force at frequency f [Hz].
k  = sqrt(2*pi*f/cf);                            % from omega = cf*k^2
k1 = 4*pi/(3*L);                                 % smallest lattice wave number
N  = ceil(2.5*k/k1) + 5;                         % modes kept (m >= n >= 0); passes 2k

% Window holds the plate plus its mirror images
x = linspace(P.x0 - 1.25*L, P.x0 + 1.25*L, P.Ng);
y = linspace(P.y0 - 1.25*L, P.y0 + 1.25*L, P.Ng);
[X, Y] = meshgrid(x, y);
inside = (Y >= 0) & (Y <= sqrt(3)*X) & (Y <= sqrt(3)*(L - X));
dA = (x(2) - x(1))*(y(2) - y(1));

W = zeros(size(X));
for m = 0:N
    for n = 0:m
        [phi, kmn] = tri_mode(m, n, k1, X, Y);
        phi0 = tri_mode(m, n, k1, P.x0, P.y0);  % <delta, phi_mn>
        if abs(phi0) < 1e-6, continue; end      % not driven by the center force

        Fd = 1;                                  % finite screw: disk average
        if kmn > 0, Fd = 2*besselj(1, kmn*P.rc)/(kmn*P.rc); end

        nrm = sum(phi(inside).^2)*dA;            % <phi_mn, phi_mn> over the plate
        A   = Fd*phi0/(nrm*(kmn^2 - (k + 1i*P.gam)^2));   % projection coefficient
        W   = W + A*phi;
    end
end

Wabs = abs(W)/max(abs(W(inside)));
Wabs(hypot(X - P.x0, Y - P.y0) < P.rc) = NaN;    % hide the screw region
end
