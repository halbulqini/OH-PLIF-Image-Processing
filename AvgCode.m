clear all;
clc;
data = dir('*.xlsx');
A=[00000];
 for i=1:length(data)
 
    filename = data(i).name; 
    B=xlsread(filename);
     B=double (B); 
     j=int2str(i);
     mx=max(max(B(1:650,:)));
% mn=min(min(c));
 y=imagesc(B,[100 mx]);
 axis image;colormap(jet);colorbar;
 title(['Image No.',j]);
 saveas(gcf,[j,'_.jpg']);
    A=double (A);
    A=A+B;
 end
 S=A./i;  % averaging
%   imshow(S,[])
%   colormap jet;
mx1=max(max(S(1:650,:)));
   filename='Avg.xlsx';
   xlswrite(filename,S);
   imagesc(S, [100 mx1]);% showing and saving the pic
axis image;colormap(jet);colorbar;
title(['AVG ']);
%  imwrite(S,'Avg.tif');
%  filename='Avg.xlsx';
%  xlswrite(filename,A);
%  imwrite(A,'Avg.tif');
 imwrite(S,jet,'Avg.jpg');
%  j=int2str(i);
%     axis image;colormap(jet);colorbar;
%     title(['Image No.',j]);
% saveas(gcf,[S,'.jpg']);
     