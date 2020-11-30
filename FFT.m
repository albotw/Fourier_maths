X = [2, 3, 4, 2];
Y = fft(X);
Z = ifft(Y);

Img = [2 1; -1 1];
Img_fft = fft2(Img);