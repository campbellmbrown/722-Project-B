function [Y_o_F] = PointReceptance(M, K, loss_factor)

force_vector = zeros(length(M), 1);
force_vector(length(force_vector) - 1) = 1;

K = K*(1 + 1i*loss_factor);

f = (0:700)*2*pi;

for i = 1:length(f)
    receptance = inv(K - (f(i)^2)*M)*force_vector;
    % Getting the second-to-last element because we are measuring at x=l,
    % and the last element is related to angle, not displacement
    Y_o_F(i) = receptance(length(receptance) - 1);
end

end