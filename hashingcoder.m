% %[Author Notes]%
% Author: Rachel Rajan
% Email : rachelrajan13@gmail.com

% %[Algorithm]%:
%This program implements the Robust Image hashing based on Ring partition and NMF algorithm presented in the follwing paper:
%Z. Tang, X. Zhang, L. Huang, Y. Dai, Robust perceptual image hashing %based on ring partition and NMF,IEEE Trans. Knowl. Data. En. 26 (3) %(2014) 711-724.
% Please reference this paper when reporting work done using this code.

%Program

clc;clear all;

%% Original Image
% take sample image peppers.png from matlab toolbox
y =imread('peppers.png');
y=rgb2ycbcr(y);
y=y(:,:,1);
y=imresize(y,[512 512]);
imshow(y);
axis on;

m=512;n=10; % n=no of rings;
x_c=m/2+0.5; % centre 
y_c=m/2+0.5;
hold on;
plot(x_c, y_c, 'g+', 'MarkerSize', 50);
r(n)=floor(m/2);% rn
A=pi*(r(n))^2;
step=A/n;
r(1)=sqrt(step/pi);%r1

%Calculating ring radius
for i=2:n-1
    r(i)=sqrt((pi*(r(i-1))^2+step)/pi);
end;
%Plotting the rings
t1 = linspace(0,2*pi,10000);
for i=1:n
     x_in = r(i)*cos(t1)+m/2+0.5;y_in = r(i)*sin(t1)+m/2+0.5;  
     hold on;
     scatter(x_in,y_in,1);
end;

% calculating Eucledian distance
for i= 1:512
    for j=1:512
    d(i,j)=sqrt((i-x_c)^2+(j-y_c)^2);
    end;
end;

op=[];
for k=1:n
    ind=find(d<r(k));% finding pixels within kth ring
    R=y(ind);%Rk vector
    u=sort(R);% uk vector
    temp=imresize(R,[step 1] );%mapping to ring area
 op=[op temp];%combining the columns
 
end;
op=double(op);
[w,h]=nnmf(op,2);
save nonmatrix h
% H=imresize(op,[512 512]);
% figure(2);
% imshow(H);

%% Rotated image

c =imread('peppers.png');
c=imrotate(c,5);
c=rgb2ycbcr(c);
c=c(:,:,1);
c=imresize(c,[512 512]);
imwrite(c,'peppers_rot1.png')
% Display the original luminance image.
figure(2)
imshow(c, []);
axis on;
m=512;
n=10; % n=no of rings;
% centre
x_c1=m/2+0.5;  
y_c1=m/2+0.5;
hold on;
plot(x_c1, y_c1, 'g+', 'MarkerSize', 50);

r1(n)=floor(m/2);% rn
A1=pi*(r1(n))^2;%Area of each ring
step1=A1/n;%Average area
r1(1)=sqrt(step1/pi);%r1

%calculating ring radius
for i=2:n-1
    r1(i)=sqrt((pi*(r1(i-1))^2+step1)/pi);
end;
%Plotting the rings
t11 = linspace(0,2*pi,10000);
for i=1:n
     x_in1 = r1(i)*cos(t11)+m/2+0.5;y_in1 = r1(i)*sin(t11)+m/2+0.5;  
     hold on;
     scatter(x_in1,y_in1,1);
end;

% calculating Eucledian distance
for i= 1:512
    for j=1:512
    d1(i,j)=sqrt((i-x_c1)^2+(j-y_c1)^2);
    end;
end;

op1=[];
for k=1:n
    ind1=find(d1<r1(k));% finding pixels within kth ring
    R1=c(ind1);%Rk vector
    u1=sort(R1);% uk vector
    temp1=imresize(R1,[step1 1] );%mapping to ring area
 op1=[op1 temp1];%combining the columns
 
end;
op1=double(op1);
[w1,h1]=nnmf(op1,2);
save nonmatrix h1

%% Tampered Image

d =imread('peppers_tamp.png');
d=rgb2ycbcr(d);
d=d(:,:,1);
imwrite(d,'peppers_tamp1.png');
% Display the original luminance image.
figure(3)
imshow(d, []);
axis on;
m=512;n=10; % n=no of rings;
% centre
x_c2=m/2+0.5;  
y_c2=m/2+0.5;
hold on;
plot(x_c2, y_c2, 'g+', 'MarkerSize', 50);

r2(n)=floor(m/2);% rn
A2=pi*(r2(n))^2;%Area of each ring
step2=A2/n;%Average area
r2(1)=sqrt(step2/pi);%r1

%calculating ring radius
for i=2:n-1
    r2(i)=sqrt((pi*(r2(i-1))^2+step2)/pi);
end;
%Plotting the rings
t12 = linspace(0,2*pi,10000);
for i=1:n
     x_in2 = r2(i)*cos(t12)+m/2+0.5;y_in2 = r2(i)*sin(t12)+m/2+0.5;  
     hold on;
     scatter(x_in2,y_in2,1);
end;

% calculating Eucledian distance
for i= 1:512
    for j=1:512
    d2(i,j)=sqrt((i-x_c2)^2+(j-y_c2)^2);
    end;
end;

op2=[];
for k=1:n
    ind2=find(d2<r2(k));% finding pixels within kth ring
    R2=d(ind2);%Rk vector
    u2=sort(R2);% uk vector
    temp2=imresize(R2,[step2 1] );%mapping to ring area
 op2=[op2 temp2];%combining the columns
 
end;
op2=double(op2);
[w2,h2]=nnmf(op2,2);
save nonmatrix h2

%% JPEG Compression

f1 = @(block_struct) dct2(block_struct.data);
f2 = @(block_struct) idct2(block_struct.data);
e = imread('peppers.png');
e=imresize(e,[512 512]);
e=rgb2gray(e);
J=blockproc(e, [8 8], f1);
depth = find(abs(J) < 20);
J(depth)= zeros(size(depth));
K = blockproc(J, [8 8], f2)/255;
imwrite(K,'peppers_jpeg.png')
Compression_ratio = numel(J)/numel(depth);

% Display the original luminance image.
figure(4)
imshow(K, []);
axis on;
m=512;n=10; % n=no of rings;
% centre
x_c3=m/2+0.5;  
y_c3=m/2+0.5;
hold on;
plot(x_c3, y_c3, 'g+', 'MarkerSize', 50);

r3(n)=floor(m/2);% rn
A3=pi*(r3(n))^2;%Area of each ring
step3=A3/n;%Average area
r3(1)=sqrt(step3/pi);%r1

%calculating ring radius
for i=2:n-1
    r3(i)=sqrt((pi*(r3(i-1))^2+step3)/pi);
end;
%Plotting the rings
t13 = linspace(0,2*pi,10000);
for i=1:n
     x_in3 = r3(i)*cos(t13)+m/2+0.5;y_in3 = r3(i)*sin(t13)+m/2+0.5;  
     hold on;
     scatter(x_in3,y_in3,1);
end;

% calculating Eucledian distance
for i= 1:512
    for j=1:512
    d3(i,j)=sqrt((i-x_c3)^2+(j-y_c3)^2);
    end;
end;

op3=[];
for k=1:n
    ind3=find(d3<r3(k));% finding pixels within kth ring
    R3=K(ind3);%Rk vector
    u3=sort(R3);% uk vector
    temp3=imresize(R3,[step3 1] );%mapping to ring area
 op3=[op3 temp3];%combining the columns
 
end;

op3=double(op3);
[w3,h3]=nnmf(op3,2);
save nonmatrix h3


%% Correaltion Coefficient Calculation

%Threshold
T=0.98;

% Hash Values(h,h1,h2,h3)

%% Analysis b/w original image and rotated image

[p,err] = polyfit(h,h1,1);   % First order polynomial
h1_fit = polyval(p,h,err);   % Values on a line
h1_dif = h1 - h1_fit;          % h1 value difference (residuals)
SSdif = sum(h1_dif.^2);      % Sum square of difference
SStot = (length(h1)-1)*var(h1);   % Sum square of h1 taken from variance
S_1 = 1-SSdif/SStot;
if S_1>T
    msgbox('Images are same')
else
    msgbox('Images are Different')
end;
    
%% Analysis b/w original image and tampered image

[p,err] = polyfit(h,h2,1);   % First order polynomial
h2_fit = polyval(p,h,err);   % Values on a line
h2_dif = h2 - h2_fit;          % h2 value difference (residuals)
SSdif = sum(h2_dif.^2);      % Sum square of difference
SStot = (length(h2)-1)*var(h2);   % Sum square of h2 taken from variance
S_2 = 1-SSdif/SStot;
if S_2>T
    msgbox('Images are same')
else
    msgbox('Images are Different')
end;

%% Analysis b/w original image and compressed image


[p,err] = polyfit(h,h3,1);   % First order polynomial
h3_fit = polyval(p,h,err);   % Values on a line
h3_dif = h3 - h3_fit;          % h3 value difference (residuals)
SSdif = sum(h3_dif.^2);      % Sum square of difference
SStot = (length(h3)-1)*var(h3);   % Sum square of h3 taken from variance
S_3 = 1-SSdif/SStot;
if S_3>T
    msgbox('Images are same')
else
    msgbox('Images are Different')
end;







