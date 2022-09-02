clear;
clc;
a=xlsread('Avg.xlsx');
% imshow(a,[]); colormap jet;
L=length(a);
b=zeros(L,L);

%Input the position of max. values at Avg. photo.
x1=124;
x2=x1;
y1=187;
y2=670;
b(x1,y1)=1801;
b(x2,y2)=1906;

%Linear gradual
d=(b(x2,y2)-b(x1,y1))/(y2-y1);
for i=y1+1:y2-1
    b(x1,i)=b(x1,y1)-d*(i-y1);
end
for i=1:y1-1
    b(x1,i)=b(x1,y1)+d*(y1-i);
end
for i=y2:L
    b(x1,i)= b(x1,y1)-d*(i-y1);
end
for i=1:L
    for j=1:L
        b1(j,i)=b(x1,i);
    end
end

%Normalization
    max1=max(max(b1));
    bn=b1./max1;
    %imshow(bn,[]); colormap jet;
  
%Dividing on the previous bpm to get the gaussian distribution
   prof=xlsread('Avgbpm1cut.xlsx');
   max2=max(max(prof));
   profn=prof./max2;
   newb=profn./bn;
   newb1=newb./max(max(newb));
   imshow(newb1,[]); colormap jet;
%Saving as Excel sheet
   filename=('bpm_mod.xlsx');
   xlswrite(filename,newb1);