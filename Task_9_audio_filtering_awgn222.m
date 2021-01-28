clear all;
close all;
clc;
[e,rs] = audioread('minions.wav');
% sound(e,rs);
% figure(1);
% plot(e);
% title('PLOT FOR ORIGINAL AUDIO');
snr = 0:1:40;
snrdb=10.^(snr/10);
%******Normalization****
[r,c]= size(e);
for i=1:c
    x=e(:,i);
    z(:,i)= x/max(x);
end
            
% figure(2);
% plot(z);
% title('PLOT FOR NORMALIZED DATA');
for j=1:length(snr)
    u=(1/snrdb(j));
    w=sqrt(u/2);
    noise=w*((randn(size(z))));
    N2=e+noise;

end
% figure(3);
% plot(N2);
% title('PLOT FOR NOISY AUDIO');
noise_audio = 'N2.wav';
audiowrite(noise_audio,N2,rs);
[N2,rs] = audioread(noise_audio);
% sound(N2,rs);
magnitude_N = abs(fft(N2));
% figure(6);
% plot(magnitude_N);

norm_N = length(magnitude_N);
figure(4);
plot([0:1/(norm_N/2-1):1], magnitude_N(1:norm_N/2));
title('PLOT FOR NORMALISED FREQUENCY');
%*******************using low pass filters*******************
[b2,a2] = butter(2,0.03,'low');
audio_filter2 = filter(b2,a2,N2);
[b,a] = ellip(1,2,30,0.03,'low');
audio_filter = filter(b,a,N2);
%******************using stopband filters********************
% [b,a] = cheby1(2,3,[0.04 0.999],'stop');
% audio_filter = filter(b,a,N2);
% [b2,a2] = cheby2(5,35,[0.04 0.999],'stop');
% audio_filter2 = filter(b2,a2,N2);
% [b,a] = ellip(1,2,30,[0.04 0.999],'stop');
% audio_filter = filter(b,a,audio_filter2);

figure(5);
plot(audio_filter);
title('PLOT FOR FILTERED AUDIO');
filter_audio = 'audio_filter.wav';
audiowrite(filter_audio,audio_filter,rs);
[audio_filter,rs] = audioread(filter_audio);
sound(audio_filter,rs);