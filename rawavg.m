clear all;
clc;

data = dir('*.tif');
C=[00000];
 for t=1:length(data)
    filename = data(t).name; 
    B=imread(filename);
    B=double(B);

C=C+B;
 end
 av=C./t;
  filename='rawavg.xlsx';
 xlswrite(filename,av);
