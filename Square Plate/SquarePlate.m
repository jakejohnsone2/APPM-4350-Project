clear
close all
clc
%% chladni_square.m -- Chladni figure on a free-edge square plate, center drive
% Solves  lap(W) + (k + i*gamma)^2 W = -delta(x-x0)delta(y-y0),  dW/dn = 0
% on 0<=x,y<=L with a Fourier cosine series (Report Sec. 2.4, Eqs. 6-7).
% Light regions in the plot = |W| ~ 0 = where the sand collects.



%% 1. Plate properties and frequency -> wave number (Sec. 2.2, Eq. 4)
% Initializing our constants of the lab plate
E = 69e9;  nu = 0.33;  rho = 2700;  h = .64e-3;   % aluminum
L = 0.1334;                                          % side length [m]
D  = E*h^3/(12*(1 - nu^2));                        % bending stiffness
cf = sqrt(D/(rho*h));                              % omega = cf*k^2  (from k^4 = rho*h*omega^2/D)

% Initializing which mode we want to see
our_n = 4;
our_m = 4;

% Creating a function to find the resonance frequency using the value of n
% and m
fnm = @(n,m) cf*(pi/L)^2*(n.^2 + m.^2)/(2*pi);     % [Hz]

% Inputing our desire n and m to find the resonance frequency of the mode
f     = fnm(our_n,our_m);              % [Hz]

% Calculating the frequency
omega = 2*pi*f;

% Calculating the k value for our system
k     = (rho*h*omega^2/D)^(1/4);
%% 2. Numerical settings (Sec. 2.7)
N    = 40;                     % modes kept per direction (n,m = 0..N)
Ng   = 500;                    % grid points per side
gam  = 0.01*pi/L;              % damping gamma [1/m]; keeps resonance finite
x0   = L/2;  y0 = L/2;         % center drive location
rc   = 0.0026;                 % screw radius [m]

%% 3. Series solution (Eq. 7)
% Creating our grid coordinates for both x and y
x = linspace(0, L, Ng);

% Creating the n and values up to our max n
n = 0:N;

% Initializing the constants
c = [1, 0.5*ones(1, N)];                         % c_0 = 1, c_n = 1/2

% Creating a meshgrid for each respective n and m pair
[NN, MM] = meshgrid(n, n);                       % columns = n, rows = m

% Calculating the value of k for each n and m pair
knm2 = (pi/L)^2*(NN.^2 + MM.^2);

% Finite screw correction
xd = sqrt(knm2)*rc;                              % finite clamp: disk-average of the source


Fd = ones(size(xd));  nz = xd > 0;               % mean-value factor 2*J1(x)/x for Helmholtz modes
Fd(nz) = 2*besselj(1, xd(nz))./xd(nz);
phi0 = cos(NN*pi*x0/L).*cos(MM*pi*y0/L).*Fd;     % phi_nm seen by the screw
wt   = L^2*(c(:)*c);                             % L^2 c_n c_m
A    = phi0 ./ (wt .* (knm2 - (k + 1i*gam)^2));  % coefficients A(m,n)
C = cos(x(:)*n*pi/L);                            % Ng x (N+1) cosine basis
W = C*A*C.';                                     % W(y_j, x_i) = sum A_mn cos cos


%% 4. Plot |W| (light = no motion = sand)
Wabs = abs(W)/max(abs(W(:)));
[Xg, Yg] = meshgrid(x);
Wabs(hypot(Xg - x0, Yg - y0) < rc) = NaN;        % hide the screw region
figure;
imagesc(x, x, Wabs, 'AlphaData', ~isnan(Wabs)); axis image xy;
colormap(flipud(gray)); caxis([0 0.15]);         % saturate so nodal lines stand out
title(sprintf('Square plate, f = %.1f Hz, Mode: n = %.0f m = %.0f', f,our_n,our_m));
xlabel('x [m]'); ylabel('y [m]');
saveas(gcf, ["square_pattern_n"+ num2str(our_n) + "_m" + num2str(our_m) + ".png"]);


%% 5. Finding Closest Mode to our Measured Value
f_meas = 3201;
T = closest_mode_square(f_meas,175,5,0);