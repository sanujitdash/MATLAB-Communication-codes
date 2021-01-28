clear all;
close all;
clc;
[e,rs] = audioread('minions.wav');
sound(e,rs);
figure(1);
plot(e);
title('PLOT FOR ORIGINAL AUDIO');
snr = 0:1:40;
snrdb=10.^(snr/10);
%******Normalization****
[r,c]= size(e);
for i=1:c
    x=e(:,i);
    z(:,i)= x/max(x);
end    
figure(2);
plot(z);
title('PLOT FOR NORMALIZED DATA');

for k=1:length(z)
    if z(k)>0
        y(k)=1;
    else
        y(k)=-1;
    end
end
for j=1:length(snr)
    u=(1/snrdb(j));
    w=sqrt(u/2);
    noise=w*(randn(1,length(y)));
    N4=y+noise;
 
end
figure(3);
plot(N4);
title('PLOT FOR NOISY AUDIO');
noise_audio = 'N4.wav';
audiowrite(noise_audio,N4,rs);
[N4,rs] = audioread(noise_audio);
% sound(N4,rs);
magnitude_N = abs(fft(N4));
norm_N = length(magnitude_N);

figure(4);
plot([0:1/(norm_N/2-1):1], magnitude_N(1:norm_N/2));
title('PLOT FOR NORMALISED FREQUENCY');
[b,a] = butter(2,0.0372,'low');   %check
audio_filter4 = filter(b,a,N4);
figure(5);
plot(audio_filter4);
title('PLOT FOR FILTERED AUDIO');
filter_audio = 'audio_filter4.wav';
audiowrite(filter_audio,audio_filter4,rs);
[audio_filter4,rs] = audioread(filter_audio);
sound(audio_filter4,rs);
