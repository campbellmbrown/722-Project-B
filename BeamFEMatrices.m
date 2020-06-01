function [K_beam,M_beam]=BeamFEMatrices(rho,area,EI,lelm)

% Matlab code to generate the element mass and stiffness matrices
% 
% DOFs = [(w, theta)_left   (w, theta)_right]
%
% *** Data required ***
% 
% rho = density
% area = cross-sectional area
% EI = bending stiffness
% lelem = length of the element
%
% ************************************

% Mass matrix "Melem"
M_beam=(rho*area*lelm/420)...
    .*[156,     22*lelm,    54,      -13*lelm;
       22*lelm, 4*lelm^2,   13*lelm, -3*lelm^2;
       54,      13*lelm,    156,     -22*lelm;
      -13*lelm, -3*lelm^2, -22*lelm,  4*lelm^2];


% Stiffness matrix "Kelem"
K_beam=(EI/(lelm^3))...
    .*[12,     6*lelm,   -12,      6*lelm;
       6*lelm, 4*lelm^2, -6*lelm,  2*lelm^2;
      -12,    -6*lelm,    12,     -6*lelm;
       6*lelm, 2*lelm^2, -6*lelm,  4*lelm^2];
