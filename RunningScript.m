clear;clc;

%% Model parameters =======================================================

L = 0.35;       % Length (m)
w = 0.02;       % Width (m)
t = 0.002;      % Thickness (m)
rho = 7850;     % Density (kg/m^3)
E = 200e9;      % Young's Modulus (pascals)
M_t = 0.02;     % Point mass (kg)

% Load theoretical mode shapes and natural frequencies
load('X');
load('Y_r');
load('omega_r');
% Given the anayltical Y = sinh(b*L)*sin(b*x) + sin(b*L)*sinh(b*x), where b
% = 0 for no harmonic excitation, the modeshape at w_nat = 0 is 0.
Y_r = [zeros(1, length(Y_r)); Y_r];
omega_r = [0, omega_r];

%% Rayleigh-Ritz method ===================================================

[nat_freqs_RR, mode_shapes_RR] = RayleighRitz(L, w, t, rho, E, M_t);
RR_xrange = linspace(0, L, length(mode_shapes_RR));
% Signs to flip mode shape
RR_sign = [1 -1 -1 -1];

for i = 1:4
    figure('Name', ['Mode ' num2str(i)])
    hold on
    plot(X, Y_r(i,:),'--');
    plot(RR_xrange, RR_sign(i)*mode_shapes_RR(i,:));
    title(['Mode shape ' num2str(i)])
    ylabel('Amplitude')
    xlabel('Distance along beam (m)')
    legend('Theoretical', 'Rayleigh-Ritz')
    ylim([-1 1]);
    grid on
end

%% Finite element method ==================================================

[nat_freqs_FE, mode_shapes_FE, M, K] = FiniteElement(L, w, t, rho, E, M_t, 200);

FE_xrange = linspace(0, L, size(mode_shapes_FE, 1));

% Signs to flip mode shape
FE_sign = [1 1 -1 1];


for i = 1:4
    figure('Name', ['Mode ' num2str(i)])
    hold on
    plot(X, Y_r(i,:),'--');
    plot(FE_xrange, FE_sign(i)*mode_shapes_FE(:,i));
    title(['Mode shape ' num2str(i)])
    ylabel('Amplitude')
    xlabel('Distance along beam (m)')
    legend('Theoretical', 'FE (200 elements)')
    ylim([-1 1]);
    grid on
end

%% Natural frequencies ====================================================
th_Hz = omega_r/2/pi;
rr_Hz = nat_freqs_RR/2/pi;
fe_Hz = nat_freqs_FE/2/pi;

fprintf('Natural Frequencies (Hz)\n');
fprintf('Theoretical:\n');
fprintf('%.2f %.2f %.2f %.2f\n', th_Hz(1),  th_Hz(2), th_Hz(3), th_Hz(4));
fprintf('Reyleigh-Ritz method:\n');
fprintf('%.2f %.2f %.2f %.2f\n', rr_Hz(1),  rr_Hz(2), rr_Hz(3), rr_Hz(4));
fprintf('Finite element:\n');
fprintf('%.2f %.2f %.2f %.2f\n', fe_Hz(1),  fe_Hz(2), fe_Hz(3), fe_Hz(4));

%% Point Receptance - via finite element method ===========================
loss_factor = 0.02;
[YoF_FE] = PointReceptance(M, K, loss_factor);
load('Y_o_F_analytical');

figure('Name','Point Receptance')
semilogy(0:700, abs(Y_o_F));    % Analytical receptance
hold on
semilogy(0:700, abs(YoF_FE));   % Finite element receptance
title('Point receptance (x = 1)')
ylabel('$\frac{Y(l)}{F}$','Interpreter', 'latex')
xlabel('Frequency (Hz)')
legend('Analytical', 'Finite element');
ylim([-inf 1])
grid on

