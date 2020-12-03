signal_origine = [0, 1, 0, 0];
signal_FFT = fft(signal_origine);
signal_iFFT = ifft(signal_FFT);

matrice_origine = [1 0; 1 0];
matrice_FFT = fft2(matrice_origine);
matrice_iFFT = ifft2(matrice_FFT);