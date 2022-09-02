clear all;
clc;
data = dir('*.tif');
A=[00000];
 for i=1:length(data)
 
    filename = data(i).name; 
    B=imread(filename);
     B=double (B);                           
     A=double (A);
    A=A+B;
 end
 S=A./i;  % averaging
  imshow(S,[])
  colormap jet;
   filename='Avgbpm1.xlsx';
   xlswrite(filename,S);
   imagesc(S);% showing and saving the pic
axis image;colormap(jet);colorbar;
title(['Avgbpm1 ']);
%  imwrite(S,'Avg.tif');
%  filename='Avg.xlsx';
%  xlswrite(filename,A);
%  imwrite(A,'Avg.tif');
%  imwrite(A,jet,'Avg.jpg');
%  j=int2str(i);
%     axis image;colormap(jet);colorbar;
%     title(['Image No.',j]);
 %saveas(gcf,[S,'.jpg']);
     