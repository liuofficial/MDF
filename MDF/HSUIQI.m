function out = HSUIQI(ref,tar)
% UIQI - calls the method described in "A Universal Image Quality Index"
% by Zhou Wang and Alan C. Bovik
n_bands = size(ref,3);
q_band = zeros(1, n_bands);
for idx1=1:n_bands
    q_band(idx1)=img_qi(ref(:,:,idx1), tar(:,:,idx1), 32);
end
out = mean(q_band);
end

function [quality, quality_map] = img_qi(img1, img2, block_size)

%========================================================================
%
%Copyright (c) 2001 The University of Texas at Austin
%All Rights Reserved.
% 
%This program is free software; you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation; either version 2 of the License, or
%(at your option) any later version.
% 
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
% 
%The GNU Public License is available in the file LICENSE, or you
%can write to the Free Software Foundation, Inc., 59 Temple Place -
%Suite 330, Boston, MA 02111-1307, USA, or you can find it on the
%World Wide Web at http://www.fsf.org.
%
%Author  : Zhou Wang 
%Version : 1.0
% 
%The authors are with the Laboratory for Image and Video Engineering
%(LIVE), Department of Electrical and Computer Engineering, The
%University of Texas at Austin, Austin, TX.
%
%Kindly report any suggestions or corrections to zwang@ece.utexas.edu
%
%Acknowledgement:
%The author would like to thank Mr. Umesh Rajashekar, the Matlab master
%in our lab, for spending his precious time and giving his kind help
%on writing this program. Without his help, this program would not
%achieve its current efficiency.
%
%========================================================================
%
%This is an efficient implementation of the algorithm for calculating
%the universal image quality index proposed by Zhou Wang and Alan C. 
%Bovik. Please refer to the paper "A Universal Image Quality Index"
%by Zhou Wang and Alan C. Bovik, published in IEEE Signal Processing
%Letters, 2001. In order to run this function, you must have Matlab's
%Image Processing Toobox.
%
%Input : an original image and a test image of the same size
%Output: (1) an overall quality index of the test image, with a value
%            range of [-1, 1].
%        (2) a quality map of the test image. The map has a smaller
%            size than the input images. The actual size is
%            img_size - BLOCK_SIZE + 1.
%
%Usage:
%
%1. Load the original and the test images into two matrices
%   (say img1 and img2)
%
%2. Run this function in one of the two ways:
%
%   % Choice 1 (suggested):
%   [qi qi_map] = img_qi(img1, img2);
%
%   % Choice 2:
%   [qi qi_map] = img_qi(img1, img2, BLOCK_SIZE);
%
%   The default BLOCK_SIZE is 8 (Choice 1). Otherwise, you can specify
%   it by yourself (Choice 2).
%
%3. See the results:
%
%   qi                    %Gives the over quality index.
%   imshow((qi_map+1)/2)  %Shows the quality map as an image.
%
%========================================================================

if (nargin == 1 | nargin > 3)
   quality = -Inf;
   quality_map = -1*ones(size(img1));
   return;
end

if (size(img1) ~= size(img2))
   quality = -Inf;
   quality_map = -1*ones(size(img1));
   return;
end

if (nargin == 2)
   block_size = 8;
end

N = block_size.^2;
sum2_filter = ones(block_size);

img1_sq   = img1.*img1;
img2_sq   = img2.*img2;
img12 = img1.*img2;

img1_sum   = filter2(sum2_filter, img1, 'valid');
img2_sum   = filter2(sum2_filter, img2, 'valid');
img1_sq_sum = filter2(sum2_filter, img1_sq, 'valid');
img2_sq_sum = filter2(sum2_filter, img2_sq, 'valid');
img12_sum = filter2(sum2_filter, img12, 'valid');

img12_sum_mul = img1_sum.*img2_sum;
img12_sq_sum_mul = img1_sum.*img1_sum + img2_sum.*img2_sum;
numerator = 4*(N*img12_sum - img12_sum_mul).*img12_sum_mul;
denominator1 = N*(img1_sq_sum + img2_sq_sum) - img12_sq_sum_mul;
denominator = denominator1.*img12_sq_sum_mul;

quality_map = ones(size(denominator));
index = (denominator1 == 0) & (img12_sq_sum_mul ~= 0);
quality_map(index) = 2*img12_sum_mul(index)./img12_sq_sum_mul(index);
index = (denominator ~= 0);
quality_map(index) = numerator(index)./denominator(index);

quality = mean2(quality_map);
end

%%%%%%%%%%%%%% Q2n aux. function
function q = onions_quality(dat1,dat2,size1)

dat1=double(dat1);
dat2=double(dat2);
dat2=cat(3,dat2(:,:,1),-dat2(:,:,2:end));
[~,~,N3]=size(dat1);
size2=size1;

% Block normalization
for i=1:N3
  [a1,s,t]=norm_blocco(squeeze(dat1(:,:,i)));
  dat1(:,:,i)=a1;
  clear a1
  if s==0
      if i==1
        dat2(:,:,i)=dat2(:,:,i)-s+1;
      else
        dat2(:,:,i)=-(-dat2(:,:,i)-s+1);   
      end
  else
      if i==1
        dat2(:,:,i)=((dat2(:,:,i)-s)/t)+1;
      else
        dat2(:,:,i)=-(((-dat2(:,:,i)-s)/t)+1);    
      end
  end
end

m1=zeros(1,N3);
m2=zeros(1,N3);

mod_q1m=0;
mod_q2m=0;
mod_q1=zeros(size1,size2);
mod_q2=zeros(size1,size2);

for i=1:N3
    m1(i)=mean2(squeeze(dat1(:,:,i)));
    m2(i)=mean2(squeeze(dat2(:,:,i)));
    mod_q1m=mod_q1m+(m1(i)^2);
    mod_q2m=mod_q2m+(m2(i)^2);
    mod_q1=mod_q1+((squeeze(dat1(:,:,i))).^2);
    mod_q2=mod_q2+((squeeze(dat2(:,:,i))).^2);
end

mod_q1m=sqrt(mod_q1m);
mod_q2m=sqrt(mod_q2m);
mod_q1=sqrt(mod_q1);
mod_q2=sqrt(mod_q2);

termine2 = (mod_q1m*mod_q2m);
termine4 = ((mod_q1m^2)+(mod_q2m^2));
int1=(size1*size2)/((size1*size2)-1)*mean2(mod_q1.^2);
int2=(size1*size2)/((size1*size2)-1)*mean2(mod_q2.^2);
termine3=int1+int2-(size1*size2)/((size1*size2)-1)*((mod_q1m^2)+(mod_q2m^2));

mean_bias=2*termine2/termine4;
if termine3==0
    q=zeros(1,1,N3);
    q(:,:,N3)=mean_bias;
else
    cbm=2/termine3;
    qu=onion_mult2D(dat1,dat2);
    qm=onion_mult(m1,m2);
    qv=zeros(1,N3);
    for i=1:N3
        qv(i)=(size1*size2)/((size1*size2)-1)*mean2(squeeze(qu(:,:,i)));
    end
    q=qv-(size1*size2)/((size1*size2)-1)*qm;
    q=q*mean_bias*cbm;
end

end
%%%%%%%%%%%%%% Q2n aux. function
function ris=onion_mult(onion1,onion2)

N=length(onion1);

if N>1
  
    L=N/2;

    a=onion1(1:L);
    b=onion1(L+1:end);
    b=[b(1),-b(2:end)];
    c=onion2(1:L);
    d=onion2(L+1:end);
    d=[d(1),-d(2:end)];


    if N==2
        ris=[a*c-d*b,a*d+c*b];
    else
        ris1=onion_mult(a,c);
        ris2=onion_mult(d,[b(1),-b(2:end)]); %%
        ris3=onion_mult([a(1),-a(2:end)],d); %%
        ris4=onion_mult(c,b);

        aux1=ris1-ris2;
        aux2=ris3+ris4;

        ris=[aux1,aux2];
    end
   
else
    ris = onion1*onion2;
end

end
%%%%%%%%%%%%%% Q2n aux. function
function ris = onion_mult2D(onion1,onion2)

[~,~,N3]=size(onion1);

if N3>1
   
    L=N3/2;

    a=onion1(:,:,1:L);
    b=onion1(:,:,L+1:end);
    b=cat(3,b(:,:,1),-b(:,:,2:end));
    c=onion2(:,:,1:L);
    d=onion2(:,:,L+1:end);
    d=cat(3,d(:,:,1),-d(:,:,2:end));


    if N3==2
        ris=cat(3,a.*c-d.*b,a.*d+c.*b); 
    else
        ris1=onion_mult2D(a,c);
        ris2=onion_mult2D(d,cat(3,b(:,:,1),-b(:,:,2:end)));
        ris3=onion_mult2D(cat(3,a(:,:,1),-a(:,:,2:end)),d);
        ris4=onion_mult2D(c,b);

        aux1=ris1-ris2;
        aux2=ris3+ris4;

        ris=cat(3,aux1,aux2);
    end
    
else
    ris = onion1.*onion2;   
end

end
%%%%%%%%%%%%%% Q2n aux. function
function [y,a,c] = norm_blocco(x)

a=mean2(x);
c=std2(x);

if(c==0)
	c = eps;
end

y=((x-a)/c)+1;

end