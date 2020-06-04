clear;

%% Model parameters =======================================================

L = 0.35;       % Length (m)
w = 0.02;       % Width (m)
t = 0.002;      % Thickness (m)
rho = 7850;     % Density (kg/m^3)
E = 200e9;      % Young's Modulus (pascals)
M_t = 0.02;     % Point mass (kg)

load('Y_r_analytical');     % Analytical mode shapes
load('omega_r_analytical'); % Analytical natural frequencies
load('X_analytical');       % X-range for analytical point receptance
load('Y_o_F_analytical');   % Analytical point receptance

%% Rayleigh-Ritz method ===================================================

[nat_freqs_RR, mode_shapes_RR] = RayleighRitz(L, w, t, rho, E, M_t);
RR_xrange = linspace(0, L, length(mode_shapes_RR));
% Signs to flip mode shape
RR_sign = [1 -1 -1 -1];

% Plotting
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
FE_sign = [-1 1 -1 1];

% Plotting
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
th_Hz = omega_r/2/pi;       % Theoretical nat. freq. (Hz)
rr_Hz = nat_freqs_RR/2/pi;  % Rayleigh-Ritz nat. freq. (Hz)
fe_Hz = nat_freqs_FE/2/pi;  % Finite element nat. freq. (Hz)

fprintf('NATURAL FREQUENCIES (Hz)\n');
fprintf('Theoretical:\n');
fprintf('%.2f %.2f %.2f %.2f\n', th_Hz(1),  th_Hz(2), th_Hz(3), th_Hz(4));
fprintf('Reyleigh-Ritz method:\n');
fprintf('%.2f %.2f %.2f %.2f\n', rr_Hz(1),  rr_Hz(2), rr_Hz(3), rr_Hz(4));
fprintf('Finite element:\n');
fprintf('%.2f %.2f %.2f %.2f\n', fe_Hz(1),  fe_Hz(2), fe_Hz(3), fe_Hz(4));

%% Point Receptance - via finite element method ===========================

loss_factor = 0.02;
f = 0:700;
[YoF_FE] = PointReceptance(M, K, f, loss_factor);

% Plotting
figure('Name','Point Receptance')
semilogy(f, abs(Y_o_F));
hold on
semilogy(f, abs(YoF_FE));
title('Point receptance (x = l)')
ylabel('$\frac{Y(l)}{F}$','Interpreter', 'latex')
xlabel('Frequency (Hz)')
legend('Analytical', 'Finite element (200 elements)');
ylim([-inf 1])
grid on