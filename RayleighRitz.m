function [natFreqs, modeShapes] = RayleighRitz(L, w, t, rho, E, m)

    A_cs = w*t;
    I = w*(t^3)/12;
    
    n = 6; % Number of segments to split the beam into
    m_s = rho*A_cs*L/n;
    F = 1;
    L = 1;
    
    for i = 1:n
       b = (2*(i - 1) + 1)*L/2/n;
       delta = F*(b^3)/(3*E*I);
       theta = F*(b^2)/2/E/I;
       for j = 1:n
          if j < i
              x = (2*(j - 1) + 1)*L/2/n;
              A(j, i) = F/6/E/I*x^2*(3*b - x);
          else
            A(j, i) = delta + (j - i)*L/n*theta;
          end
       end
    end
    natFreqs = 1;
    modeShapes = 1;
end