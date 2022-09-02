% P1=reading the first phoro.
% E is to put each photo inside it before summing.
% M is a zeros matrix to sum all photos in it.
% we use double precision to make it available to have large values of cells and get correct values of average.
% AV is the average one single shot of flat flame.
% S is the matrix used to make Beam Profile Matrix.
% m is to exclude the high wrong values at the bottom of AV.
% F is a matrix to put the coloumns of S in it and 'fill' the other values.
clc;
clear;

P1=imread('B00001.tif');
L=length(P1);
E=zeros(L,L);
data = dir('*.tif')
M=P1.*0;
M=double(M);
for i=1:200
    filename = data(i).name; 
    B=imread(filename);
    E=double(B);
    M=M+E;
end
AV=M./i;
xlswrite('Avg1.xlsx',AV);
% xfile=xlsread('Avg1.xlsx');
% imwrite(xfile,jet,'Avgt.jpg');

z=1;
for i=150:50:650
    for j=1:L
        S(j,z)=AV(j,i);
    end
    z=z+1;
end
S=double(S);
m=S(1:740,:);
maax=max(m);

 % Smoothing with 5 pixels
 for i=1:11
  N(1:L,i)=S(1:L,i)./maax(i);
 end
 BeamProfile= N;
 for x=1:11
 smoothed(1,x)=BeamProfile(1,x);
 smoothed(2,x)=BeamProfile(2,x);
 smoothed(767,x)=BeamProfile(767,x);
 smoothed(768,x)=BeamProfile(768,x);
 for i=3:1:766
 smoothed(i,x)=(BeamProfile(i-2,x)+BeamProfile(i-1,x)+BeamProfile(i,x)+BeamProfile(i+1,x)+BeamProfile(i+2,x))/5;
 end
 end
 
 % Smoothing with 3 pixels
%   for x=1:11
%  smoothed(1,x)=BeamProfile(1,x);
%  smoothed(768,x)=BeamProfile(768,x);
%  for i=2:1:767
%  smoothed(i,x)=(BeamProfile(i-1,x)+BeamProfile(i,x)+BeamProfile(i+1,x))/3;
%  end
%  end

 % Smoothing with 7 pixels
%  for i=1:11
%   N(1:L,i)=S(1:L,i)./max(i);
%  end
%  BeamProfile= N;
%  for x=1:11
%  smoothed(1,x)=BeamProfile(1,x);
%  smoothed(2,x)=BeamProfile(2,x);
%  smoothed(3,x)=BeamProfile(3,x);
%  smoothed(766,x)=BeamProfile(766,x);
%  smoothed(767,x)=BeamProfile(767,x);
%  smoothed(768,x)=BeamProfile(768,x);
%  for i=4:1:765
%  smoothed(i,x)=(BeamProfile(i-3,x)+BeamProfile(i-2,x)+BeamProfile(i-1,x)+BeamProfile(i,x)+BeamProfile(i+1,x)+BeamProfile(i+2,x)+ BeamProfile(i+3,x))/7;
%  end
%  end

%  imshow(smoothed,[]);
%  colormap jet;


     % Filling with dx between 9 coloumns
%  F=zeros(L,L);
%  F=double(F);
%  % Gradual filling from the beginning 
% %  for i=1:150
% %      for j=1:L
% %      dx=smoothed(j,1)/150;
% %      F(j,i)=dx*i;
% %      end
% %  end
% 
% % The same value from 1 to 150
%     for i=1:150
%          F(:,i) = smoothed(:,1);
%     end
% 
%      for i=151:200
%          for j=1:L
%      dx=(smoothed(j,2)-smoothed(j,1))/50;
%      F(j,i)=dx*(i-150)+F(j,150);
%          end
%      end
%      for i=201:250
%          for j=1:L
%             dx=(smoothed(j,3)-smoothed(j,2))/50;
%             F(j,i)=dx*(i-200)+F(j,200);
%          end
%      end
%      for i=251:300
%          for j=1:L
%             dx=(smoothed(j,4)-smoothed(j,3))/50;
%             F(j,i)=dx*(i-250)+F(j,250);
%          end
%      end
%      for i=301:350
%          for j=1:L
%             dx=(smoothed(j,5)-smoothed(j,4))/50;
%             F(j,i)=dx*(i-300)+F(j,300);
%          end
%      end
%      for i=351:400
%          for j=1:L
%             dx=(smoothed(j,6)-smoothed(j,5))/50;
%             F(j,i)=dx*(i-350)+F(j,350);
%          end
%      end
%      for i=401:450
%          for j=1:L
%             dx=(smoothed(j,7)-smoothed(j,6))/50;
%             F(j,i)=dx*(i-400)+F(j,400);
%          end
%      end
%      for i=451:500
%          for j=1:L
%             dx=(smoothed(j,8)-smoothed(j,7))/50;
%             F(j,i)=dx*(i-450)+F(j,450);
%          end
%      end
%      for i=501:550
%          for j=1:L
%             dx=(smoothed(j,9)-smoothed(j,8))/50;
%             F(j,i)=dx*(i-500)+F(j,500);
%          end
%      end
%      for i=551:600
%          for j=1:L
%             dx=(smoothed(j,10)-smoothed(j,9))/50;
%             F(j,i)=dx*(i-550)+F(j,550);
%          end
%      end
%      for i=601:650
%          for j=1:L
%             dx=(smoothed(j,11)-smoothed(j,10))/50;
%             F(j,i)=dx*(i-600)+F(j,600);
%          end
%      end
%      
%      % Gradual filling to the end
% %      for i=651:L
% %          for j=1:L
% %          dx=smoothed(j,11)/(L-650);
% %          F(j,i)=-dx*(i-650)+F(j,650);
% %          end
% %      end
% 
% % The same value from 650 to the end
%     for i=651:L
%          F(:,i) = smoothed(:,11);
%     end
%     
%      imshow(F,[]);
%      colormap jet;
%      xlswrite('9&FillF.xlsx',F);



     % Filling with arith. mean of 9 coloumns
%      NewS=smoothed(:,2:10);
%      oneS=sum(NewS,2)/9;
%      for i=1:L
%          NewF(:,i) = oneS;
%      end
%      imshow(NewF,[]);
%      colormap jet;
%      xlswrite('9F.xlsx',NewF);



   % Filling with arith. mean of 500 coloumns 
  allS = AV(:,150:650);
  m2=allS(1:740,:);
    maax2=max(m2);
 for i=1:500
  N2(1:L,i)=allS(1:L,i)./maax2(i);
 end
 BeamProfile2= N2;
 
%  imshow(BeamProfile2,[]);
%  colormap jet;
 
  for x=1:500
 smoothed2(1,x)=BeamProfile2(1,x);
 smoothed2(2,x)=BeamProfile2(2,x);
 smoothed2(767,x)=BeamProfile2(767,x);
 smoothed2(768,x)=BeamProfile2(768,x);
 for i=3:1:766
 smoothed2(i,x)=(BeamProfile2(i-2,x)+BeamProfile2(i-1,x)+BeamProfile2(i,x)+BeamProfile2(i+1,x)+BeamProfile2(i+2,x))/5;
 end
  end
 
     oneallS=sum(smoothed2,2)/500;
     for i=1:L
         NewallF(:,i) = oneallS;
     end
     imshow(NewallF,[]);
     colormap jet;
     
     xlswrite('500F.xlsx',NewallF);
    
   xlswrite(filename,NewallF);
 imwrite(NewallF,'500F.tif');
 imwrite(NewallF,jet,'500F.jpg');