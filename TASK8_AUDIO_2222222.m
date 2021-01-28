clear all;
close all;
clc;
[e,rs] = audioread('minions.wav');
% sound(e,rs);
% figure(1);
% plot(e);
% title('PLOT FOR ORIGINAL AUDIO');
k = 25000;
%******Normalization****
[r,c]= size(e);
for i=1:c
    x=e(:,i);
    z(:,i)= x/max(x);
end
            
% figure(2);
% plot(z);
% title('PLOT FOR NORMALIZED DATA');

% encrypt_audio = floor(z*k);
encrypt_audio = floor(z*k);
degrade_audio ='encrypt_audio.wav';
% figure(3);
% plot(encrypt_audio);
% title('PLOT FOR ENCRYPTED AUDIO');
audiowrite(degrade_audio,encrypt_audio,rs);
[encrypt_audio,rs] = audioread(degrade_audio);
% sound(encrypt_audio,rs);

decrypt_audio = encrypt_audio/k;
received_audio = 'decrypt_audio.wav';
% figure(4);
% plot(decrypt_audio);
% title('PLOT FOR DECRYPTED AUDIO');
audiowrite(received_audio,decrypt_audio,rs);
[decrypt_audio,rs] = audioread(received_audio);
% sound(decrypt_audio,rs);

