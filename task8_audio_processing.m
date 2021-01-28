clear all;
close all;
clc;
[e,rs] = audioread('CantinaBand3.wav');
% sound(e,rs);
figure(1);
plot(e);
title('PLOT FOR ORIGINAL AUDIO');

k = 25000;
%******Normalization****
[r,c]=size(e);
for i=1:c   % for i=1:size(e,2)
    x=e(:,i);
    z(:,i)= x/max(x);
end
% max_val = max(e);
% min_val = abs(min(e));
% if max_val(1,1)>min_val(1,1)
%     max_n(1,1)=max_val(1,1);
% else
%     max_n(1,1)=min_val(1,1);
% end
% % if max_val(1,2)>min_val(1,2)
% %     max_n(1,2)=max_val(1,2);
% % else
% %     max_n(1,2)=min_val(1,2);
% % end
% norm_ip(:,1) = e(:,1)./max_n(1,1);
% norm_ip(:,2) = e(:,2)./max_n(1,2);
figure(2);
plot(z);
title('PLOT FOR NORMALIZED DATA');

en_audio = floor(z*k);
degrade_audio ='en_audio.wav';
figure(3);
plot(en_audio);
title('PLOT FOR ENCRYPTED AUDIO');
% audiowrite(degrade_audio,en_audio,rs);
% [en_audio,rs] = audioread(degrade_audio);
% sound(en_audio,rs);

de_audio = en_audio/k;
received_audio = 'de_audio.wav';
figure(4);
plot(de_audio);
title('PLOT FOR DECRYPTED AUDIO');
audiowrite(received_audio,de_audio,rs);
[de_audio,rs] = audioread(received_audio);
sound(de_audio,rs);


