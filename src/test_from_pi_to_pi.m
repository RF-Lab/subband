clc, clf, clear all;

curPath = pwd() ;
cd('.\tsim') ;
modelPath = pwd() ;
cd( curPath ) ;
addpath(modelPath) ;

A = 3; E = A^2 / 2;
fs = 5;
fd = 20 ;
N = 20;

SNR_dB = -27:1:10 ;
sigma = E ./ (10 .^ (SNR_dB./10)) ;
SNR = E ./ sigma ;

%w1 = -pi; w2 = pi;
freq_range = 0.1; w1 = (fs - freq_range)/fd * 2 *pi; w2 = (fs + freq_range)/fd * 2 *pi ;
freq_range_2 = 1; w1_2 = (fs + freq_range_2)/fd * 2 *pi; w2_2 = (fs + 2*freq_range_2)/fd * 2 *pi ;
%omega_range=0.5*pi; w1_2 = (fs)/fd * 2 *pi  - omega_range; w2_2 = (fs)/fd * 2 *pi + 3*omega_range;

% get sub band coefficients
F0 = get_sb_matrix_1(N, w1, w2, 0) 
F1 = get_sb_matrix_1(N, w1, w2, 1) 
F2 = get_sb_matrix_1(N, w1, w2, 2) 

F0_2 = get_sb_matrix_1(N, w1_2, w2_2, 0) 
F1_2 = get_sb_matrix_1(N, w1_2, w2_2, 1) 
F2_2 = get_sb_matrix_1(N, w1_2, w2_2, 2) 

F0_pi = get_sb_matrix_1(N, -pi, pi, 0) 
F1_pi = get_sb_matrix_1(N, -pi, pi, 1) 
F2_pi = get_sb_matrix_1(N, -pi, pi, 2) 

experiment_size = 1000;
freq_acf = zeros(length(sigma), 1);
var_acf = zeros(length(sigma), 1);
freq_sub = zeros(length(sigma), 1);
var_sub = zeros(length(sigma), 1);
freq_sub_2 = zeros(length(sigma), 1);
var_sub_2 = zeros(length(sigma), 1);
freq_sub_pi = zeros(length(sigma), 1);
var_sub_pi = zeros(length(sigma), 1);

for snr_range = 1:length(sigma)
    for num=1:experiment_size
        s = E * cos(2*pi*fs/fd * (0:N-1));
        x = s + sqrt(sigma(snr_range)) * randn(1, length(s));

        %% acf
        X = fft(x, 1*N);
        XX = X .* conj(X);
        acf_full = ifft(XX);
        %acf_full = [x.'*x, x.' * circshift(x, 1), x.' * circshift(x, 2)] / N;

        %plot(acf_full(1:20));
        acf = ar_model([acf_full(1); acf_full(2); acf_full(3)]) ;
        [poles1, omega0_acf, Hjw0_1] = get_ar_pole(acf) ;
        fs_acf_est = omega0_acf*fd/2/pi;
        freq_acf(snr_range) = freq_acf(snr_range) + (fs_acf_est)^2;
        var_acf(snr_range) = var_acf(snr_range) + (fs_acf_est - fs)^2;

        % sub band
        x = x.';
        r = [x.' * F0 * x, x.' * F1 * x, x.' * F2 * x] / (2*pi);

        sub = ar_model([r(1); r(2); r(3)]) ;
        [poles1, omega0_sub, Hjw0_1] = get_ar_pole(sub) ;
        fs_sub_est = omega0_sub*fd/2/pi;
        freq_sub(snr_range) = freq_sub(snr_range) + (fs_sub_est)^2;
        var_sub(snr_range) = var_sub(snr_range) + (fs_sub_est - fs)^2;

        % sub band 2
        r_2 = [x.' * F0_2 * x, x.' * F1_2 * x, x.' * F2_2 * x] / (2*pi);

        sub_2 = ar_model([r_2(1); r_2(2); r_2(3)]) ;
        [poles1, omega0_sub_2, Hjw0_1] = get_ar_pole(sub_2) ;
        fs_sub_est = omega0_sub_2*fd/2/pi;
        freq_sub_2(snr_range) = freq_sub_2(snr_range) + (fs_sub_est)^2;
        var_sub_2(snr_range) = var_sub_2(snr_range) + (fs_sub_est - fs)^2;
        
        % sub from -pi to pi
        r_pi = [x.' * F0_pi * x, x.' * F1_pi * x, x.' * F2_pi * x] / (2*pi);
                
        sub_pi = ar_model([r_pi(1); r_pi(2); r_pi(3)]) ;
        [poles1, omega0_sub_pi, Hjw0_1_pi] = get_ar_pole(sub_pi) ;
        fs_sub_pi_est = omega0_sub_pi*fd/2/pi;
        freq_sub_pi(snr_range) = freq_sub_pi(snr_range) + (fs_sub_pi_est)^2;
        var_sub_pi(snr_range) = var_sub_pi(snr_range) + (fs_sub_pi_est - fs)^2;
        
    end
    
    freq_acf(snr_range)  = sqrt(freq_acf(snr_range)  / experiment_size);
    freq_sub(snr_range)  = sqrt(freq_sub(snr_range)  / experiment_size);
    freq_sub_2(snr_range)  = sqrt(freq_sub_2(snr_range)  / experiment_size);
    freq_sub_pi(snr_range)  = sqrt(freq_sub_pi(snr_range)  / experiment_size);
    
    var_acf(snr_range)  = var_acf(snr_range)  / experiment_size;
    var_sub(snr_range)  = var_sub(snr_range)  / experiment_size;
    var_sub_2(snr_range)  = var_sub_2(snr_range)  / experiment_size;
    var_sub_pi(snr_range)  = var_sub_pi(snr_range)  / experiment_size;
end

subplot(2, 1, 1),
    plot(SNR_dB, freq_acf, '-rx', SNR_dB, freq_sub, '-go',  SNR_dB, freq_sub_2, '-b+', SNR_dB, freq_sub_pi, '-kd');
    legend('acf', sprintf('sub band from %.2f Hz to %.2f Hz', w1*fd/2/pi, w2*fd/2/pi), sprintf('sub band from %.2f Hz to %.2f Hz', w1_2*fd/2/pi, w2_2*fd/2/pi),'sub band - \pi to \pi')
    xlabel('SNR, dB'), ylabel('Freq, Hz');

title(sprintf('Fs=%.2f Hz\t Fd=%.2f Hz', fs, fd));
    
subplot(2, 1, 2);
    plot(SNR_dB, var_acf, '-rx', SNR_dB, var_sub, '-go',  SNR_dB, var_sub_2, '-b+', SNR_dB, var_sub_pi, '-kd');
    legend('acf', sprintf('sub band from %.2f Hz to %.2f Hz', w1*fd/2/pi, w2*fd/2/pi), sprintf('sub band from %.2f Hz to %.2f Hz', w1_2*fd/2/pi, w2_2*fd/2/pi),'sub band - \pi to \pi')
    xlabel('SNR, dB'), ylabel('Variance');

%title(sprintf('From %.2f Hz to %.2f Hz', w1*fd/2/pi, w2*fd/2/pi))


%fprintf('Sub band:%.2f\t pure ACF:%.2f\n', freq_sub, freq_acf)

% remove model path
rmpath(modelPath) ;