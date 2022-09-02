% Correction Image for Beam Profile
% Reading Images from tiff and forming one assembled image
clear all;
clc;
data = dir('*.tif');
% for i=1:length(data)   % Reading in different formats
% filename = data(i).name;
% h=imread(filename);
% D= uint8(h);
%  j = int2str (i);
% imwrite(D,[j,'.tif']);
% end
A=[00000];
 for i=1:length(data)
 
    filename = data(i).name; 
    B=imread(filename);
     B=double (B);                           
     A=double (A);
    A=A+B;
 end
 A=A./i;  % averaging
  imshow(A,[])
  colormap jet;
 z=1;
  for i=200:1:600;     % specifying limits
       for j=1:length(A);
 c(1:j,z)=A(1:j,i);
       end
     z=z+1;
  end
  % plotting of vectors for checking  tilting 
 for i=1:401 %mod from 401 to 1001
     hold on;
  plot(c(:,i),length(A))
  end
 c=double(c);
 m=c(1:600,:);
 max1=max(m);
  for i=1:401 %mod from 401 to 1001
   N(1:length(A),i)=c(1:length(A),i)./max1(i);     % Normalization
  end
  BeamProfile= N;
  for x=1:401 %mod from 401 to 1001
  smoothed(1,x)=BeamProfile(1,x);
  smoothed(2,x)=BeamProfile(2,x);
  smoothed(767,x)=BeamProfile(767,x);
 smoothed(768,x)=BeamProfile(768,x);
 for i=3:1:766
 smoothed(i,x)=(BeamProfile(i-2,x)+BeamProfile(i-1,x)+BeamProfile(i,x)+BeamProfile(i+1,x)+BeamProfile(i+2,x))/5;
 end
 end
  imshow(smoothed,[]);
 colormap jet;
 f=[00000];
 for i=1:401 %mod from 401 to 1001
    t=smoothed(1:length(A),i);
     f=double (f);                           
     t=double (t);
    f=t+f;
 end
f=f./401;%mod from 401 to 1001
for i=1:length(A)
 F(1:length(A),i)=f;
 end
imshow(F,[]);
 colormap jet;
 filename='Correction3.xlsx';
 xlswrite(filename,F);
 imwrite(F,'correction3.tif');
 imwrite(F,jet,'correction3.jpg');