L = 0.35;
w = 0.02;
t = 0.002;
rho = 7850;
E = 200e9;
m = 0.02;
n = 4;

[natFreqs, modeShapes] = RayleighRitz(L, w, t, rho, E, m);

%[natFreqs, modeShapes] = FiniteElement();

% Load theoretical mode shapes and natural frequencies
load('X');
load('Y_r');
load('omega_r');

figure('Name', 'Mode Shapes')
for i = 1:n
    subplot(2, 2, i);
    plot(linspace(0, L, n), modeShapes(:,i));
    hold on
    plot(X, Y_r(i,:),'--');
    title(['Mode shape ' num2str(i)])
    ylabel('Amplitude')
    xlabel('Distance along beam (m)')
    grid on
end

figure('Name', 'Natural Frequencies')
plot(omega_r(1:n)/2/pi, '-or');
hold on
plot(natFreqs/2/pi, '-ok');
title('Natural Frequencies')
ylabel('Natural Frequency (Hz)')
xlabel('Mode')
xlim([0 n+1]);
grid on 