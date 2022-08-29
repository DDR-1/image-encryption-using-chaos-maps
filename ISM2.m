img=imread('lena.png');
img1 = imread('encrypted.png');
img2= imread('diff_encrypted.png');
img3 = imread('encrypted1.png');
display(['Entropy:  ',num2str(entropy(img1))])
plot_hist(img,img1);
NPCR_and_UACI(img1,img2);
corr(img);
corr(img1);
plot_diff(img1,img2);

function plot_hist(image,encr_image)
figure
subplot(1,2,1);
imhist(image);
title('Original Image');
subplot(1,2,2);
imhist(encr_image);
title('Encrypted Image');
end

function NPCR_and_UACI(i1,i2)
D = 0;
UACI = 0;
[M,N] = size(i1);

for i = 1:M
  for j = 1:N
    if i1(i,j) ~= i2(i,j)
      D = D + 1;
      UACI = UACI + abs(((double(i1(i,j))-double(i2(i,j))))/(M*N*255));
    end
   end
end    

NPCR = (D/(M*N))*100;
UACI = UACI*100;
display(['NPCR:  ',num2str(NPCR)]);
display(['UACI:  ',num2str(UACI)]);
end

function corr(img)
%Horizontal Correlation
I=double(img);
c_diag = corrcoef(I(1:end-1, 1:end-1), I(2:end, 2:end));
 c_vert = corrcoef(I(1:end-1, :), I(2:end, :));
 c_horz = corrcoef(I(:,1:end-1), I(:, 2:end));
 display(c_horz);
 display(c_vert);
 display(c_diag);
end

function plot_diff(i1,i2)
[M,N] = size(i1);
i3=zeros(M,N);
for i = 1:M
  for j = 1:N
      i3(i,j) = i1(i,j) - i2(i,j);
  end
end
figure
imshow(i3);
imwrite(i3,'diffx.png','BitDepth',8);
end  

