clear all;
close all;
clc;

[e,r] = audioread('audio1.mp3');
% name = audio1.wav;
% sound(d,r);
file_name = fopen('audio1.mp3');
d = fread(file_name, [1,Inf],'uint8');
%%d = fread(file_name',ubit1');
fclose(file_name);
b = uint8(rem(floor((1./[128, 64, 32,16,8,4,2,1].').*d),2));
b = b(:);
[r c] = size(b);
n = reshape(b,c,r);


% audio_bin = de2bi(e);
% [r c] = size(audio_bin);
% n = reshape(audio_bin,r,c);
                    
data=n;
% nn=length(n);
z=[];

snr = 0:1:20;
snrdb=10.^(snr/10);

for k=1:length(n)
    if data(k)>0
        y(k)=1;
    else
        y(k)=-1;
    end
end

for j=1:length(snr)
    for l=1:20
   % N = awgn(complex(y),snr(j));
    u=(1/snrdb(j));
    w=sqrt(u/2);
    noise=w*(randn(1,length(y)));
    N=y+noise;

        for i=1:length(n)
            if N(i)>0
                z(i) = 1;
            else
                z(i) = 0;
            end
        end
% t = reshape(z,r1,c1);
% image_dec = bi2de(t);
% image_op = reshape(image_dec,r,c);
% image_out=uint8(image_op);
% 
% 
error = sum(z~=data);
BER(l)=error/length(n);
    end
BER_avg(j)= sum(BER)/20;
    

end
BER_Theoretical=0.5*erfc(sqrt(snrdb));
figure(2);
% imshow(image_out);
% title('received image');
% figure(3);
semilogy(snr,BER_avg,'b-');
hold on;
semilogy(snr,BER_Theoretical,'r*');
grid on;
xlabel('SNR in dB--->');
ylabel('BER--->');
axis([0 12 10^-10 1]);
title ('Simulating BPSK in AWGN channel using audio as data bits');
legend ('BER Practical', 'BER Theoretical');