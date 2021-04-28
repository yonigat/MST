function [sum_sigma] = sum_sigma_calc(E0)
    kevtoj = 1.60218e-16;
    E0 = E0*kevtoj;
    re = 2.818e-13; %[cm]
    re2 = re^2;
    c = 3e8; %[m/s]
    c2 = c^2;
    me = 9.11e-31; %[kg]
    mec2 = me*c2; %[jouls]
    Z = 1;

    n_energy = 1000;
    E1min = E0*mec2/(mec2 + 2*E0);
    E1 = linspace(E1min, E0, n_energy);

    epsilon = E1/E0;

    cost = 1 - mec2./epsilon/E0 + mec2/E0;
    cos2t = cost.^2;

    sin2t = 1 - cos2t;
    
    f = 1./epsilon + epsilon;
    g = 1 - (epsilon .* sin2t) ./ (1 + epsilon.^2);
    
    dsigma = pi * re2 * (mec2 / E0) * Z .* f .* g; % need to divide by 2*pi

    sum_sigma = sum(dsigma)*((E1(2)-E1(1))/E0);
end 


