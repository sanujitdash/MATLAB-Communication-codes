  clear all;
close all;
clc;
                    
b = imread('eyes.jpg');
% figure(5);
% image(b);
% title('Original Image');
if size(b,3)==3
    b=rgb2gray(b);
end
image_grey = b;
[r c] = size(image_grey);
m = reshape (image_grey,1,numel(image_grey));  
image_bin = de2bi(m);
[r1 c1] = size(image_bin);
n = reshape(image_bin,1,numel(image_bin));
data=n;
nn=length(n);
z=[];
w=[];
v=[];
mse_qpsk = [];
mse_bpsk = [];
snr = 0:1:40;
snrdb=10.^(snr/10);
error1 = zeros(1,length(snr));
error2 = zeros(1,length(snr));
error3 = zeros(1,length(snr));

for i=1:length(n) %or use i=1:nn
    if data(i)>0
        y(i)=1;
    else
        y(i)=-1;
    end
    R(i) = y(i)+1j*y(i);
end
for k =1:length(snr)
noise = (1/sqrt(2))*(randn(1,length(y))+1j*randn(1,length(y)));
rayleigh_h = (1/sqrt(2))*(randn(1,length(y))+1j*randn(1,length(y)));
awgn_N = noise*(10^(-snr(k)/20));

%Rayleigh channel  
rayl = rayleigh_h.*R;
rayl_ch = awgn_N+rayl;  
 %Equalization   
N = rayl_ch./rayleigh_h;


rayl2 = rayleigh_h.*y;
rayl_ch2 = awgn_N+rayl2;   
N2 = rayl_ch2./rayleigh_h;

%error counting
w = real(N)>=0;
z = imag(N)>=0;
v = real(N2)>0;

t = reshape(w,r1,c1);
image_dec = bi2de(t);

t2 = reshape(z,r1,c1);%%%
image_dec2 = bi2de(t2);%%%

image_avg = (image_dec+image_dec2)/2;%%%
image_op = reshape(image_avg,r,c);
image_out = uint8 (image_op);

t3 = reshape(v,r1,c1);
image_dec3 = bi2de(t3);
image_op3 = reshape(image_dec3,r,c);
image_out3 = uint8 (image_op3);

error1(k) = sum(w~=data);
error2(k) = sum(z~=data);
error3(k) = sum(v~=data);
% Mean squared error
mse_qpsk(k) = immse(image_out,image_grey);
disp(mse_qpsk);
mse_bpsk(k) = immse(image_out3,image_grey);
disp(mse_bpsk);
end

ber1 = error1./nn;
ber2 = error2./nn;
BER_qpsk1 = (ber1+ber2)/2;
BER_bpsk1 = error3/nn;

BER_theoretical_qpsk = 0.5*(1-sqrt(snrdb./(snrdb+1)));

% figure(1);
% imshow(image_out);
% title('Received Image(QPSK)');
% figure(2);
% imshow(image_out3);
% title('Received Image(BPSK)');
figure(6);
plot(mse_qpsk,'b-');
hold on;
plot(mse_bpsk,'r*');
grid on;
xlabel('SNR (in dB)--->');
ylabel('Mean Squared Error--->');
axis([0 30 0 2100]);
title ('MEAN SQUARED ERROR');
legend ('MSE(QPSK)', 'MSE(BPSK)');
% 
% figure(3);
% semilogy(snr,BER_qpsk1,'b-','LineWidth',4);
% hold on;
% semilogy(snr,BER_theoretical_qpsk,'r*','LineWidth',4);
% grid on;
% xlabel('SNR in dB--->');
% ylabel('BER--->');
% axis([0 30 10^-7 10^0]);
% title ('Simulating BER of QPSK in Rayleigh channel using image as data bits');
% legend ('BER Simulated(QPSK)', 'BER Theoretical(QPSK)');
% 
% figure(4);
% semilogy(snr,BER_qpsk1,'g<','LineWidth',4);
% hold on;
% semilogy(snr,BER_theoretical_qpsk,'mx','LineWidth',4);
% hold on;
% semilogy (snr,BER_bpsk1,'b-','LineWidth',4);
% axis ([0 30 10^-6 10^0]);
% grid on;
% xlabel ('SNR in dB--->');
% ylabel ('BER---->');
% title ('Simulating BER of QPSK in Rayeligh Fading channel and the comparative plot using image as data bits');
% legend ('BER of QPSK Simulated(Rayleigh)','BER of QPSK Theoretical(Rayleigh)','BER of BPSK simulated(Rayleigh)','BER of BPSK Theoretical(Rayleigh)');
