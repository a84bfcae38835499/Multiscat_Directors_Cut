function I = Gaussian2D(x, y, mu, s)
    I = (1/(s*sqrt(2*pi)))*exp(-(x - mu(1)).^2/(2*s^2)).*exp(-(y - mu(2)).^2/(2*s^2));
end