clc
clear all

%Reading image
c0 = imread('lena.png');
c0 = squeeze(c0(:,:,3));

%flag and offset are for security analysis
offset1=0.000000;   %offset 1 is to change key for key sensitivity test
flag=0;     %flag is to change one random pixel for differential analysis

% Getting Image Dimensions
[M,N] = size(c0);
cd = double(c0);

% extracting all bit one by one
% since values are in form of 8-bit integer
% from c1 to c8 respectively
c1 = mod(cd, 2);
c2 = mod(floor(cd/2), 2);
c3 = mod(floor(cd/4), 2);
c4 = mod(floor(cd/8), 2);
c5 = mod(floor(cd/16), 2);
c6 = mod(floor(cd/32), 2);
c7 = mod(floor(cd/64), 2);
c8 = mod(floor(cd/128), 2);

% Creating A1 and A2 matrices containing the bitplane
A1 = [c1 c2 c3 c4 ];
A2 = [c5 c6 c7 c8 ];
L = 4*M*N;
X = zeros(1,M*N);

%Changing one bit in image for differential attack
if flag == 1
    temp = find(A1 == 1);
    r1 = randi(size(temp,1));
    A1(temp(r1)) = ~A1(temp(r1));
end

% Key for diffusion phase
x0 = 0.3+offset1;
n0 = 0.2;
% Key for confusion phase
y0 = 0.2+offset1;
u0 = 0.4;

X(1) = x0;
for i= 1:M*N-1
    X(i+1)=pwlcm(X(i),n0);
end
X = mod(floor(X*10^14),256);

x1 = mod(X, 2);
x2 = mod(floor(X/2), 2);
x3 = mod(floor(X/4), 2);
x4 = mod(floor(X/8), 2);
x5 = mod(floor(X/16), 2);
x6 = mod(floor(X/32), 2);
x7 = mod(floor(X/64), 2);
x8 = mod(floor(X/128), 2);

b1 = [x1 x2 x3 x4];
b2 = [x5 x6 x7 x8];
% initiating the encryption
% Diffusion phase
% sum calculates the sum of each column then calling sum to calculate the overall
% calling sum twice because A2 is a matrix
sum1 = sum(sum(A2));

% circle shift the matrix to the right by the value in sum1
A11 = circshift(A1,sum1(1),2);
% Initializing B1
B1 = zeros(L,1);
B1(1) = xor(xor(xor(A11(1),A11(L)),A2(1)),b1(1));

for i = 2:L
   B1(i) = xor(xor(xor(A11(i),A11(i-1)),A2(i)),b1(i));
end

% Calling sum on B1 only once since B1 is an array.
sum2 = sum(B1);

A22 = circshift(A2,sum2,2);
%A22 = [1,1,0,0,0];
% Initializing B2
B2 = zeros(L,1);
B2(1) = xor(xor(xor(A22(1),A22(L)),B1(1)),b2(1));

for i = 2:L
  B2(i) = xor(xor(xor(A22(i),A22(i-1)),B1(i)),b2(i));
end

% confusion phase

sum3 = sum2 + sum(B2); 
s0 = mod(y0+sum3/L,1);
S = zeros(1,2*L);
S(1) = s0;

for i=1:2*L-1
    if S(i)>=0 && S(i)<u0
        S(i+1)=S(i)/u0;
    end
    if S(i)>=u0 && S(i)<0.5
        S(i+1)=(S(i)-u0)/(0.5-u0);
    end
    if S(i)>=0.5 && S(i)<1-u0
        S(i+1)=(1-u0-S(i))/(0.5-u0);
    end
    if S(i)>=1-u0 && S(i)<1
        S(i+1)=(1-S(i))/u0;
    end    
end

s1 = S(1:L);
s2 = S(L+1:2*L);

Y = mod(floor(s1 * 10^14),L) + 1;
Z = mod(floor(s2 * 10^14),L) + 1;

for i = 1:L
  [B2(Y(i)),B1(i)] = deal(B1(i),B2(Y(i)));
end

for j = 1:L
  [B1(Z(j)),B2(j)] = deal(B2(j),B1(Z(j)));
end


% Simplicity sake we can assume n=1. However, to implement n rounds,
% we would have to loop the diffusion and confusion phase.

% converting the array C1  and C2 into matrices
E1 = reshape(B1,[M,4*N]);
E2 = reshape(B2,[M,4*N]);

% splitting the matrices into bitplanes
E11 = mat2cell(E1,M,[N N N N]);
E22 = mat2cell(E2,M,[N N N N]);

% extracting bitplances from E11 and E22
e1 = cell2mat(E11(1));
e2 = cell2mat(E11(2));
e3 = cell2mat(E11(3));
e4 = cell2mat(E11(4));
e5 = cell2mat(E22(1));
e6 = cell2mat(E22(2));
e7 = cell2mat(E22(3));
e8 = cell2mat(E22(4));

% combining bitplanes again to form the encrypted image
ee = (2 * (2 * (2 * (2 * (2 * (2 * (2 * e8 + e7) + e6) + e5) + e4) + e3) + e2) + e1);

% plotting original image in first subplot
subplot(1, 3, 1);
imshow(c0);
title('Original Image');

% plotting recombined image in second subplot
subplot(1,3,2);
ee=uint8(ee);
imshow(ee);
title('Encrypted Image');
if flag==0
    imwrite(ee,'encrypted.png','BitDepth',8);
end
if flag==1
    imwrite(ee,'diff_encrypted.png','BitDepth',8);
end
if offset1 ~= 0
    imwrite(ee,'encrypted1.png','BitDepth',8);
end
%DECRYPTION
d = ee;
dd = double(d);

%Extracting All Bits From the Images Into 8 Variables
d1 = mod(dd, 2);
d2 = mod(floor(dd/2), 2);
d3 = mod(floor(dd/4), 2);
d4 = mod(floor(dd/8), 2);
d5 = mod(floor(dd/16), 2);
d6 = mod(floor(dd/32), 2);
d7 = mod(floor(dd/64), 2);
d8 = mod(floor(dd/128), 2);

%Creating Matrices With Bit-Planes
D1 = [d1 d2 d3 d4 ];
D2 = [d5 d6 d7 d8 ];

%CONFUSION PHASE

%Key for the Confusion Phase
y0 = 0.2;
u0 = 0.4;

%sum3 = sum2 + sum(B2); 
s0 = mod(y0+sum3/L,1);
S = zeros(1,2*L);
S(1) = s0;

for i=1:2*L-1
    if S(i)>=0 && S(i)<u0
        S(i+1)=S(i)/u0;
    end
    if S(i)>=u0 && S(i)<0.5
        S(i+1)=(S(i)-u0)/(0.5-u0);
    end
    if S(i)>=0.5 && S(i)<1-u0
        S(i+1)=(1-u0-S(i))/(0.5-u0);
    end
    if S(i)>=1-u0 && S(i)<1
        S(i+1)=(1-S(i))/u0;
    end    
end

s1 = S(1:L);
s2 = S(L+1:2*L);

Y = mod(floor(s1 * 10^14),L) + 1;
Z = mod(floor(s2 * 10^14),L) + 1;

for j = L:-1:1
  [D2(j),D1(Z(j))] = deal(D1(Z(j)),D2(j));
end

for i = L:-1:1
  [D1(i),D2(Y(i))] = deal(D2(Y(i)),D1(i));
end

%DIFFUSION PHASE

%Key for the Diffusion Phase
x0 = 0.3;
n0 = 0.2;

X(1) = x0;
for i= 1:M*N-1
    X(i+1)=pwlcm(X(i),n0);
end
X = mod(floor(X*10^14),256);

x1 = mod(X, 2);
x2 = mod(floor(X/2), 2);
x3 = mod(floor(X/4), 2);
x4 = mod(floor(X/8), 2);
x5 = mod(floor(X/16), 2);
x6 = mod(floor(X/32), 2);
x7 = mod(floor(X/64), 2);
x8 = mod(floor(X/128), 2);

b1 = [x1 x2 x3 x4];
b2 = [x5 x6 x7 x8];

%Initializing Matrices for Solving the System of Equations
Mat1 = speye(L-1) + circshift(speye(L-1),-1,2);
Mat1(1,L-1) = 0;
Mat1(1,L) = 1;

Mat2 = zeros(L-1,1);
Mat2(1,1) = bitxor(bitxor(D1(1),D2(1),'uint8'),b2(1),'uint8');
for i=1:L-1
  Mat2(i,1) = bitxor(bitxor(D1(i),D2(i),'uint8'),b2(i),'uint8');
end

%Solving and Reshaping the Values to Form the Image Matrix
X = mod(Mat1\Mat2,2);

X22 = reshape(X,[M,4*N]);
X2 = circshift(X22,-sum2,2);

if log2(M)<5 && rem(M,2)==0 && rem(M,3)~=0
    if rem(log2(M),2)~=0 || rem(log2(M),4)==0
        X2(:) = ~X2;
    end
end

Mat3 = zeros(L-1,1);
Mat3(1,1) = bitxor(bitxor(D1(1),X2(1),'uint8'),b1(1),'uint8');
for i=1:L-1
    Mat3(i,1) = bitxor(bitxor(D1(i),X2(i),'uint8'),b1(i),'uint8');
end

Y2 = mod(Mat1\Mat3,2);

X11 = reshape(Y2,[M,4*N]);
X1 = circshift(X11,-sum1(1),2);

if rem(M,2)==0
    if floor(log2(M))~=log2(M) || log2(M)>5
        if rem(log2(M),2)~=0 || rem(log2(M),4)==0
            X1(:) = ~X1;
        end
    end
end

%Splitting the Matrices to Bit-Planes
I11 = mat2cell(X1,M,[N N N N]);
I22 = mat2cell(X2,M,[N N N N]);

%Extracting the Bit-Planes
i1 = cell2mat(I11(1));
i2 = cell2mat(I11(2));
i3 = cell2mat(I11(3));
i4 = cell2mat(I11(4));
i5 = cell2mat(I22(1));
i6 = cell2mat(I22(2));
i7 = cell2mat(I22(3));
i8 = cell2mat(I22(4));

%Combining the Bit-Planes to Form the Decrypted Image
ii = (2 * (2 * (2 * (2 * (2 * (2 * (2 * i8 + i7) + i6) + i5) + i4) + i3) + i2) + i1);

% Plotting the Decrypted Image
subplot(1,3,3);
ii=uint8(ii);
ii=imcomplement(ii);
imshow(ii);
title('Decrypted Image');



function x = pwlcm(x,n)
% x = (0,1) and n = (0,0.5)
if 0.5 < x && x <= 1
		x = pwlcm(1-x,n);
end
if 0 <= x && x < n
		x = x/n;
end
if n <= x && x <= 0.5
		x = (x-n)/(0.5-n);
end
end


