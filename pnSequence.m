close all;
clc;
flip_flop = input ('Enter the number of flip flops');
pn =(2^flip_flop)-1;
a = ones(1,flip_flop);

%*****Alternate method for calculating 1's and 0's*****
% one_c=0;
% zero_c=0;
%******************************************************
for i=1:pn
    z=a;
    p(i)=0;
    p(i) = xor(z(1,(flip_flop-1)),z(1,flip_flop));
    a(1,1) = p(i);
    for j=1:(flip_flop-1)
        a(1,(j+1))=z(1,j);
    end
%*****Alternate method for calculating 1's and 0's*****
%     if p(i)==1
%         one_c=one_c+1;
%     else
%         zero_c=zero_c+1;
%     end
%***************************************************
end
disp(p);
%*******************BALANCE PROPERTY****************
one_c=(sum(p==1));   %Uncomment the above counting process if not
zero_c=(sum(p==0));  %interested in using this one and comment these.
if one_c==zero_c+1
   disp ('Balance Property Satisfied');
else
   disp ('Balance Property NOT Satisfied');
end
%*************AUTOCORRELATION PROPERTY*************
for i=1:pn
    if p(i)==0
        p1(i)=-1;
    else
        p1(i)=1;
    end
end
s2=[];
for k=-1:pn+1
    s=circshift(p1,k);
    s1=sum(p1.*s)/pn;
    s2=[s2 s1];
end
figure;
k=-1:pn+1;
plot(k,s2,'b-');
xlabel ('Time Lag---->');
ylabel ('Autocorrelation---->');
axis ([-2 9 -0.2 1]);
legend ('Calculation of Autocorrelation');
