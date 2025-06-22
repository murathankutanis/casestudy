%%%%%%%%%%%%%%%%%%%%% Murathan Kutanis - Teknik Çalışma 1 %%%%%%%%%%%%%%%%%

% 1) 3 GHz CW (0 dBm)
% 2) %60 dB SNR ile WGN ekleme
% 3) –40 dBc (2. harm.), –50 dBc (3. harm.) ekleme
% 4) Güç spektrumu
% 5) 750 MHz merkezli 0 dBm LFM darbe
% 6) PRF=10 kHz, Pw=40 us → duty cycle
% 7) Sweep: fc–10→fc+9 MHz, 1 MHz adım, 20 frekans
% 8) Mixer: CW × adımlı LFM
% 9) 3750 MHz merkezli, BW=100 MHz, derecesi ≤5 band-pass filtre
% -------------------------------------------------------------------------

%% Params
Z0    = 50;            % Ohm
P_dBm =   0;           % dBm
P     = 10^(P_dBm/10)/1e3;   % W
A_cw  = sqrt(2*P*Z0);       % Amplitude (V)
fc1   =  3e9;          % CW frequeny
Fs    = 20e9;          % Sampling Frequency (≥2·max(fc))
Ts    = 1/Fs;          
N     = 4096;          % Spectrum Resolution
t_cw  = (0:N-1)*Ts;    % Time Vector

%% 1) CW
s1 = A_cw*cos(2*pi*fc1*t_cw);

%% 2) White Gaussian Noise
P_sig  = mean(s1.^2); %Average power
SNRdB  = 60;
P_noise= P_sig/10^(SNRdB/10);
n      = sqrt(P_noise)*randn(size(s1));
s2     = s1 + n;

%% 3) 2. ve 3. harmonik 
A2 = A_cw*10^(-40/20);    % 2
A3 = A_cw*10^(-50/20);    % 3
s3 = s2 + A2*cos(2*pi*2*fc1*t_cw) + A3*cos(2*pi*3*fc1*t_cw);

%% 4) Güç spektrumu
S3 = fftshift(fft(s3));
f  = (-N/2:N/2-1)*(Fs/N);
Pwr= abs(S3).^2/N;
figure(1);
plot(f/1e9, 10*log10(Pwr));
xlabel('Frekans [GHz]'); ylabel('Güç [dB]');
title('Güç Spektrumu');
grid on

%% 5) LFM 
fc2     = 750e6;            % Merkez frekans [Hz]
A_pulse = A_cw;             % Tepe genlik [V]
tau     = 40e-6;            % Darbe genişliği [s]
BW      = 20e6;             % Örnek bant genişliği [Hz]
t_p     = 0:Ts:tau-Ts;      % Darbe zaman ekseni

% Chirp oranı:
K       = BW / tau;         % Hz/s

phi     = 2*pi*(fc2*t_p + 0.5*K*t_p.^2);
s5      = A_pulse * cos(phi);

figure;
plot(t_p*1e6, s5);
xlabel("Zaman [µs]");
ylabel("Genlik [V]");
title(sprintf("LFM Darbe"));
grid on;

%% 6) PRF ve duty cycle
PRF   = 10e3;            % Hz
Trep  = 1/PRF;           % s
duty  = tau/Trep;        
fprintf('Duty cycle = %.1f%%\n', duty*100);

%% 7) Sweep
f_start    = fc2 - 10e6;
f_end      = fc2 +  9e6;
freq_steps = linspace(f_start, f_end, 20);

% Adımlı LFM sinyali: her adım için tek bir cw pulse
s7 = zeros(1, round(20*Trep*Fs));
for k=1:20
    idx0 = round((k-1)*Trep*Fs)+1;
    idx1 = idx0 + numel(t_p)-1;
    s7(idx0:idx1) = A_pulse*cos(2*pi*freq_steps(k)*t_p);
end

t_tr = (0:length(s7)-1) * Ts;

figure;
plot(t_tr*1e3, s7);
xlabel('Time [ms]');
ylabel('Amplitude [V]');
title('Step‐Frequency CW Pulse Train (20 steps)');
grid on;

%% 8) Mixer
L = min(numel(s3),numel(s7));
s8 = s3(1:L) .* s7(1:L);

t8 = (0:L-1) * Ts; 
figure;
plot(t8*1e3, s8);            % ms cinsinden gösterim
xlabel('Time [ms]');
ylabel('Amplitude [V]');
title('Mixer Çıkışı: s8 = s3 \times s7');
grid on;

%% 9) Band-pass filtre
fc_filt   = 3750e6;      % merkez frekans
BW        = 100e6;
f_low     = fc_filt - BW/2;
f_high    = fc_filt + BW/2;
ord       = 5;           

% normalize edilmiş kesme frekansları [0,1]
Wn = [f_low f_high]/(Fs/2);

b  = fir1(ord, Wn, 'bandpass');
s9 = filter(b,1,s8);

S9 = fftshift(fft(s9, N));
P9 = abs(S9).^2/N;
figure(2);
plot(f/1e9, 10*log10(P9));
xlabel('Frekans [GHz]'); ylabel('Güç [dB]');
title('Filtrelenmiş Radar İşareti Güç Spektrumu');
grid on
