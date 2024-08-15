clear; close all; clc;
%% Plot things! wowie!

data = readmatrix("diffrac10001.out",NumHeaderLines=7,FileType="text");
diffracChannels = int32(data(:,2:3));
diffracValues = data(:,4);
%disp("min(diffractionChannels(:,1),[],'all') = ")
%disp(min(diffractionChannels(:,1),[],"all"))
%class((min(diffractionChannels(:,1),[],"all")))
minChX = (min(diffracChannels(:,1),[],"all"));
maxChX = (max(diffracChannels(:,1),[],"all"));
minChY = (min(diffracChannels(:,2),[],"all"));
maxChY = (max(diffracChannels(:,2),[],"all"));

nChX = maxChX - minChX;
nChY = maxChY - minChY;

plotValues = NaN(size(diffracValues,1),1);
plotCoordsX = NaN(size(diffracValues,1),1);
plotCoordsY = NaN(size(diffracValues,1),1);

%we want something like
%[channelX,channelY,val] -> [x,y,val] -> plot!

%Now we need to make the reciprocal lattice and plot them ig or some shit
%idk. since we know the diffraction channels, I suppose it's a simple case
%of:
% - multiply the diffraction channel by the lattice vector
% - make two matricies storing x and y of those
% - plot the data against that! ez!

a1=[const.c*sqrt(3)/2,-const.c/2];
a2=[0,const.c];
a3=[0,0,const.c];
%A1 = a1;
%A2 = a2;
[b1,b2,b3] = Reciprocal([a1,0],[a2,0],a3);
[b1,b2] = Reciprocal(a1,a2);

X = NaN(nChX,nChY);
Y = NaN(nChX,nChY);

for chx = minChX:maxChX
    for chy = minChY:maxChY
      i = chx-minChX+1;
      j = chy-minChY+1;
      %disp("chx = " + chx + ", chy = " + chy);
      %disp("i   = " + i + ", j   = " + j);
      %disp("X = " + string(B1(1)*double(chx) + B2(1)*double(chy)))
      %disp("Y = " + string(B1(2)*double(chx) + B2(2)*double(chy)))
      X(i,j) = B1(1)*double(chx) + B2(1)*double(chy);
      Y(i,j) = B1(2)*double(chx) + B2(2)*double(chy);
    end
end

indx = 1;

for chx = minChX:maxChX
    for chy = minChY:maxChY
      i = chx-minChX+1;
      j = chy-minChY+1;
      %disp("chx = " + chx + ", chy = " + chy);
      %disp("i   = " + i + ", j   = " + j);
      %disp("X = " + string(B1(1)*double(chx) + B2(1)*double(chy)))
      %disp("Y = " + string(B1(2)*double(chx) + B2(2)*double(chy)))
      for k = 1:size(diffracChannels,1)
        %disp("diff= " + diffracChannels(k,1) + ", " + diffracChannels(k,2))
        %disp("chan= " + chx + ", " + chy)
        result = diffracChannels(k,:) == [chx,chy];
        if(result)
          disp("Found channel!!! chan= " + chx + ", " + chy)
          plotValues(indx) = diffracValues(k);
          plotCoordsX(indx) = X(i,j);
          plotCoordsY(indx) = Y(i,j);
          indx = indx + 1;
        else
          continue
        end
      end
    end
end
scatter(plotCoordsX,plotCoordsY,'filled','o',SizeData=max(plotValues,1.e-15)*5000)
daspect([1,1,1])