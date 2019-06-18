function Wdem = dwhann(dem)

    [Ny, Nx] = size(dem);

    % detrend and or background subtract
    dem = dem - min(min(dem));
    dem(isnan(dem)) = 0;

    % Hann-windowed DEM
    a=(Nx-1)/2; b=(Ny-1)/2;
    W = zeros(size(dem));

    for m=1:Nx
        for n=1:Ny
            phi = atan( (n-b) / (m-a) );
            r = sqrt( (m-a)^2+(n-b)^2 );
            r_prime = sqrt( (a^2)*(b^2) / ((b^2) * (cos(phi)^2) + (a^2)*(sin(phi)^2)));
            if r<=r_prime
                W(n,m) = 0.5*(1+cos(pi*r/r_prime));
            else
                W(n,m) = 0;
            end
        end
    end
    
    Wdem = W .* dem;
 
end