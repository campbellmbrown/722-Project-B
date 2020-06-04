function [YoF_FE] = PointReceptance(M, K, f, loss_factor)
% PointReceptance function
% Returns the point receptance, excited at x = L and measured at x = L
% INPUTS ======================
% M - finite element mass matrix
% K - finite element stiffness matrix
% f - frequency range (Hz)
% loss_factor - loss factor of beam
% OUTPUTS =====================
% Y0F_FE - finite element point receptance
% AUTHOR ======================
% Campbell Brown - 738509729

force_vector = zeros(length(M), 1);
% Force is applied at the end of the beam (corresponding to the
% second-to-last DOF)
force_vector(length(force_vector) - 1) = 1;

% Making K complex by including loss factor
K = K*(1 + 1i*loss_factor);

% Converting frequencies to rad/s
f_rad = f*2*pi;

for i = 1:length(f_rad)
    receptance = inv(K - (f_rad(i)^2)*M)*force_vector;
    % Getting the second-to-last element because we are measuring at x = L,
    % and the second-to-last element corresponds to displacement DOF
    YoF_FE(i) = receptance(length(receptance) - 1);
end

end