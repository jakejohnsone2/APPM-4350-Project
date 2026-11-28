%% SquarePlate.m -- Chladni figures on a free-edge square plate, center drive
% Model (Report Sec. 2.3-2.4):
%   lap(W) + (k + i*gamma)^2 W = -delta(x-x0) delta(y-y0),   dW/dn = 0 on the edges
% solved with the Fourier cosine series, Eq. (10).
%
% The center screw is modeled as a point force (averaged over the screw head).
% A mode is only driven if it is nonzero at the screw, phi_nm(x0,y0) ~= 0, which at
% the center means n and m both even. Every pattern below is computed at a MEASURED
% lab frequency and labeled with the nearest mode the center force can drive.
%
% Light regions in the plot = |W| ~ 0 = where the sand collects.
% Needs closest_mode_square.m in the same folder.
clear; clc; close all

%% 1. Plate properties and frequency -> wave number (Sec. 2.2)
% Setting/Calculating our initial variables
E   = 69e9;  nu = 0.33;  rho = 2700;  h = 0.64e-3;   % aluminum
L   = 0.1334;                                        % side length [m]
D   = E*h^3/(12*(1 - nu^2));                         % bending stiffness
cf  = sqrt(D/(rho*h));                               % omega = cf*k^2

%% 2. Numerical settings (Sec. 2.7
% Setting numerical settings so simulation isn't too big
P.N   = 40;                 % modes kept per direction (n,m = 0..N)
P.Ng  = 500;                % grid points per side
P.gam = 0.01*pi/L;          % damping gamma [1/m]; keeps resonance finite
P.x0  = L/2;  P.y0 = L/2;   % center drive
P.rc  = 0.0026;             % screw head radius [m]

%% 3. Measured frequencies from the lab [Hz]
% Initializing our frequencies that we found when doing the lab
f_lab = [175 1608 2565 3201 4020];

%% 4. Pattern at each measured frequency
% For loop that goes through each of our frequencies within the lab
for f = f_lab
    % Finding the response given the respective frequency
    [W, x] = square_response(f, L, cf, P);
    % Finding the closest mode to our response
    M = closest_mode_square(f, L, cf, P.x0, P.y0, 1);   % nearest driven mode

    % Taking the magnitude of our complex field
    Wabs = abs(W)/max(abs(W(:)));
    % Creating the coordinate grid
    [Xg, Yg] = meshgrid(x);
    Wabs(hypot(Xg - P.x0, Yg - P.y0) < P.rc) = NaN;      % hide the screw region

    % Plotting for the respective frequency
    figure;
    imagesc(x, x, Wabs, 'AlphaData', ~isnan(Wabs)); axis image xy;
    % Making our zero amplitudes white so the nodal response stands out
    colormap(flipud(gray)); caxis([0 0.15]);
    title(sprintf('Square, f = %.0f Hz (nearest driven mode (%d,%d), %.0f Hz)', ...
          f, M.n, M.m, M.f_pred_Hz));
    xlabel('x [m]'); ylabel('y [m]');
    saveas(gcf, sprintf('square_pattern_f%dHz.png', round(f)));
end

%% 5. Final Table
% Creating the table containing the modes, lab freq, and predicted freq
T = closest_mode_square(f_lab, L, cf, P.x0, P.y0, 1);
disp(T);

%% -------------------------------Functions--------------------------------
function [W, x] = square_response(f, L, cf, P)
%SQUARE_RESPONSE  |W| on an Ng x Ng grid for a center force at frequency f [Hz].
    k = sqrt(2*pi*f/cf);                             % from omega = cf*k^2
    if P.N*pi/L < 2*k
        warning('N = %d may be too small for k = %.1f 1/m. Increase N.', P.N, k);
    end
    % Initializing the x coordinate scheme
    x = linspace(0, L, P.Ng);

    % Setting the number of sums to our max n value
    n = 0:P.N;

    % Setting our coefficients
    c = [1, 0.5*ones(1, P.N)];                       % c_0 = 1, c_n = 1/2

    % Creating a meshgrid of the different m and n
    [NN, MM] = meshgrid(n, n);                       % columns = n, rows = m

    % Calculating k for our different ns and ms
    knm2 = (pi/L)^2*(NN.^2 + MM.^2);
    
    % Finite screw: average the point force over a disk of radius rc.
    % For a Helmholtz mode this multiplies phi(x0,y0) by 2*J1(k*rc)/(k*rc).
    xd = sqrt(knm2)*P.rc;
    Fd = ones(size(xd));  nz = xd > 0; 
    Fd(nz) = 2*besselj(1, xd(nz))./xd(nz);
    
    % Calculating our coefficients
    phi0 = cos(NN*pi*P.x0/L).*cos(MM*pi*P.y0/L).*Fd;   % <delta, phi_nm>
    nrm  = L^2*(c(:)*c);                               % <phi_nm, phi_nm> = L^2 c_n c_m
    A    = phi0 ./ (nrm .* (knm2 - (k + 1i*P.gam)^2)); % projection coefficients, Eq. (9)
    
    % Calculating the eigenfunctions
    C = cos(x(:)*n*pi/L);                            % Ng x (N+1) cosine basis

    % Calculating the response
    W = C*A*C.';                                     % W(y_j, x_i) = sum A_mn cos cos
end
