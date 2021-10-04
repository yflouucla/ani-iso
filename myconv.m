function h = myconv(f,g)


[ny,nx]=size(f);
H = psf2otf(g, [ny,nx]);

h = ifft2(fft2(f).*H);

h = real(h);
