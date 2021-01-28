clear all;
close all;
clc;
N_sub = 128;  % No of subcarriers
cp_len = 1;  % Cyclic prefix length
Ts = 1e-3;    % Sampling period of channel
df = 0;       % Max Doppler frequency shift
np = 4;       % No of pilot symbols
m = 4;         % No of symbols for QAM modulation
ofdm_frames = 10^3; % No of OFDM frames
data = round((m-1)*rand((N_sub-2*np),ofdm_frames));
const = qammod([0:m-1],m);
qam_data = qammod(data,m);
pilot_data = [zeros(np,ofdm_frames); qam_data ; zeros(np,ofdm_frames)];  % Pilot Insertion
%OFDM symbol
ifft_data = (128/sqrt(120))*ifft(pilot_data,N_sub);
Cy_pre = [ifft_data((128-cp_len+1):128,:); ifft_data];   % Cyclic prefix
[r,c] = size(Cy_pre);
Tx_Data = Cy_pre;
% Frequency selective channel with 4 taps
tau = [0 1e-5 3.5e-5 12e-5];          % Path delays
pdb = [0 -1 -1 -3];                   % Avg path power gains
h = rayleighchan(Ts, df, tau, pdb);
h.StoreHistory = 0;
h.StorePathGains = 1;
h.ResetBeforeFiltering = 1;
% SNR of channel
EbNo = 0:1:40;
EsNo= EbNo + 10*log10(120/128)+ 10*log10(128/144);      % symbol to noise ratio
snr= EsNo - 10*log10(128/144); 
% Transmit through channel
ber = zeros(1,length(snr));
demod_data = zeros((N_sub-2*np),ofdm_frames);
for i = 1:length(snr)
    for j = 1:c     % Transmit frame by frame
        hx = filter(h,Tx_Data(:,j).');   % Pass through Rayleigh channel
        a = h.PathGains;
        AM = h.channelFilter.alphaMatrix;
        g = a*AM;          % Channel coefficients
        G(j,:) = fft(g,N_sub);  % DFT of channel coefficients
        awgn_noise = awgn(hx,snr(i));   % Add AWGN noise
% Receiver
        Rx = awgn_noise(cp_len+1:r);                                % Removal of cyclic prefix 
        fft_data = (sqrt(120)/128)*fft(Rx,N_sub)./G(j,:);   % Frequency Domain Equalization
        demod_data(:,j) = qamdemod(fft_data(5:124),m);     % Removal of pilot and Demodulation 
    end
    ber(i) = sum(sum(demod_data~=data))/((N_sub-2*np)*ofdm_frames);
end
ber_theoretical = berfading(snr,'QAM',m,1);
figure(1);
semilogy(EbNo,ber_theoretical,'r*','LineWidth',4);
hold on;
semilogy(EbNo,ber,'b-','linewidth',4);
grid on;
xlabel('EbNo in dB--->');
ylabel('BER--->');
axis([0 20 10^-5 10^0]);
title ('Simulating BER of QAM in RAYLEIGH channel using OFDM');
legend ('BER Theoretical(QAM)','BER Simulated(QAM)');