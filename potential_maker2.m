clear; close all; clc;
rng default;
rng("shuffle");
%Number of grid points, number of Z points
Nxy = 32; Nz = 100; Nsuper = 1;
zMax = 6; zMin = -2;%units Å
fileprefix = "mos2_test_"

plotPot = true;

%a1=[const.a,0];
%a2=[0,const.a];
%a3=[0,0,const.a];

a1=[-const.c,0];
a2=[const.c/2,const.c*sqrt(3)/2];
a3=[0,0,const.c];
[b1,b2,b3] = Reciprocal([a1,0],[a2,0],a3);

%% Defect density calculations
cellArea = const.c^2 * sqrt(1-(dot(a1,a2)/(const.c^2))^2);
disp("Unit cell area = " + cellArea + "Å^2")
%% Now create pristine potentials
V = zeros(Nxy,Nxy,Nz);
Vinterp = zeros(Nxy,Nxy,Nz);
X = zeros(Nxy,Nxy);
Y = zeros(Nxy,Nxy);
Xsuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
Ysuper = zeros(Nxy*Nsuper,Nxy*Nsuper);
Vsuper = zeros(Nsuper*Nxy,Nsuper*Nxy,Nz);
Z = linspace(zMin,zMax,Nz);

for i = 0:Nxy-1
    for j = 0:Nxy-1
        X(i+1,j+1) = (a1(1)*i+a2(1)*j)./Nxy;
        Y(i+1,j+1) = (a1(2)*i+a2(2)*j)./Nxy;
    end
end

for i = 0:Nxy*Nsuper-1
    for j = 0:Nxy*Nsuper-1
        Xsuper(i+1,j+1) = (a1(1)*i+a2(1)*j)./Nxy;
        Ysuper(i+1,j+1) = (a1(2)*i+a2(2)*j)./Nxy;
    end
end

for k = 1:Nz
      %V(:,:,k) = Vfunc_LiF_Wolken(X,Y,Z(k));
      V(:,:,k) = Vfunc_MoS2(X,Y,Z(k));
end
for z = 1:Nz
    for nx = 1:Nxy:Nsuper*Nxy
        for ny = 1:Nxy:Nsuper*Nxy
            Vsuper(nx:nx+Nxy-1,ny:ny+Nxy-1,z) = V(:,:,z);
        end
    end
end

%% Import Min's DFT

importfile("DFT_Pure.mat")
x1=[const.d,0];
x2=[const.d/2,const.d * sqrt(3)/2];
DFTsuper = 2;
XDFTsuper = zeros(12*DFTsuper);
YDFTsuper = zeros(12*DFTsuper);
N1DFT = linspace(0,const.c,12*DFTsuper);
N2DFT = linspace(0,const.c,12*DFTsuper);
N1DFT_a = [N1DFT const.c*(12*DFTsuper+1/(12*DFTsuper))];
N2DFT_a = [N2DFT const.c*(12*DFTsuper+1/(12*DFTsuper))];
ZDFT = linspace(1.5,6,19);

for i = 0:12*DFTsuper-1
    for j = 0:12*DFTsuper-1
        XDFTsuper(i+1,j+1) = (x1(1)*i+x2(1)*j)./12;
        YDFTsuper(i+1,j+1) = (x1(2)*i+x2(2)*j)./12;
    end
end

XDFTsuper = XDFTsuper - const.c/(sqrt(3)); %makes the 0,0 point be a sulphur
XDFTsuper = XDFTsuper - x1(1) - x2(1);%shifts potential to overlap with the smaller unit cell, for interpolation purposes
YDFTsuper = YDFTsuper - x1(2) - x2(2);


theta = 30;
rotMat = [cosd(theta) -sind(theta);
          sind(theta)  cosd(theta)];

for i = 1:size(XDFTsuper,1)
  for j = 1:size(YDFTsuper,1)
    vIn = [XDFTsuper(i,j); YDFTsuper(i,j)];
    vOut = rotMat * vIn;
    XDFTsuper(i,j) = vOut(1);
    YDFTsuper(i,j) = vOut(2);
  end
end

VDFT = zeros(12,12,19);
for k = 1:19
  VDFT(:,:,k) = pagetranspose(Pot_M(:,:,k))*1000;
end
VDFT_a = zeros(12+1,12+1,19);
VDFT_a(1:12,1:12,:) = VDFT;
VDFT_a(12+1,1:12,:) = VDFT(1,:,:);
VDFT_a(1:12,12+1,:) = VDFT(:,1,:);
VDFT_a(end ,end ,:) = VDFT(end,end,:);

VDFTsuper = zeros(DFTsuper*12,DFTsuper*12,19);
for z = 1:19
    for nx = 1:12:DFTsuper*12
        for ny = 1:12:DFTsuper*12
            VDFTsuper(nx:nx+12-1,ny:ny+12-1,z) = pagetranspose(Pot_M(:,:,z))*1000;
        end
    end
end

VDFTsuper_a = zeros(DFTsuper*12+1,DFTsuper*12+1,19);
VDFTsuper_a(1:DFTsuper*12,1:DFTsuper*12  ,:) = VDFTsuper;
VDFTsuper_a(DFTsuper*12+1,1:DFTsuper*12  ,:) = VDFTsuper(1,:,:);
VDFTsuper_a(1:DFTsuper*12,DFTsuper*12+1  ,:) = VDFTsuper(:,1,:);
VDFTsuper_a(end          ,end            ,:) = VDFTsuper(end,end,:);

[y1,y2,y3] = Reciprocal([x1,0],[x2,0],a3);



%===
%% Now interpolate the DFT data into a useful basis
interpolateDFTdata = false;

if(interpolateDFTdata)
    VDFTvect = zeros(DFTsuper*DFTsuper*12*12*19,1);
    XDFTvect = VDFTvect;
    YDFTvect = VDFTvect;
    ZDFTvect = VDFTvect;
    index = 0;
    for k = 1:19
      z = ZDFT(k);
       for j = 1:12*DFTsuper
        for i = 1:12*DFTsuper
          if(index + 1 ~= 144*DFTsuper*DFTsuper*(k-1)+12*DFTsuper*(j-1)+i)
            error("F")
          end
          index = 144*DFTsuper*DFTsuper*(k-1)+12*DFTsuper*(j-1)+i;
          %disp("index = " + num2str(index))
          XDFTvect(index) = XDFTsuper(i,j);
          YDFTvect(index) = YDFTsuper(i,j);
          ZDFTvect(index) = z;
          VDFTvect(index) = VDFTsuper(i,j,k);
        end
      end
    end
    InterpolatedFn = scatteredInterpolant(XDFTvect,YDFTvect,ZDFTvect,VDFTvect,'natural','none');
    Xvect = squeeze(zeros(Nxy*Nxy*Nz,1));
    Yvect = Xvect;
    Zvect = Xvect;
    for k = 1:Nz
      z = Z(k);
       for j = 1:Nxy
        for i = 1:Nxy
          index2 = Nxy*Nxy*(k-1)+Nxy*(j-1)+i;
          %disp("index2 = " + num2str(index2))
          Xvect(index2) = X(i,j);
          Yvect(index2) = Y(i,j);
          Zvect(index2) = z;
        end
      end
    end
    %Vvect = interp3(XDFTvect,YDFTvect,ZDFTvect,VDFTsuper,Xvect,Yvect,Zvect,'linear');
    Vvect = InterpolatedFn(Xvect,Yvect,Zvect);%<- Beware! this step takes absolutely forever
    if(anynan(Vvect))
      error("Nan found!")
    end  
    for k = 1:Nz
      for j = 1:Nxy
        for i = 1:Nxy
          Vinterp(i,j,k) = Vvect(Nxy*Nxy*(k-1)+Nxy*(j-1)+i);
        end
      end
    end
  Vinterpsuper = zeros(Nsuper*Nxy,Nsuper*Nxy,Nz);
  for z = 1:Nz
      for nx = 1:Nxy:Nsuper*Nxy
          for ny = 1:Nxy:Nsuper*Nxy
              Vinterpsuper(nx:nx+Nxy-1,ny:ny+Nxy-1,z) = Vinterp(:,:,z);
          end
      end
  end
end

%writematrix(Vsuper,"V.csv")

%Vsuper = readmatrix("V_boyao.csv");
%Vsuper = reshape(Vsuper,[Ncell,Ncell,Nz]);

%===
%%
copyInterp = false;
if(copyInterp)
  Vsuper = Vinterpsuper;
end

potStructArray(1).V = Vsuper;
potStructArray(1).a1=Nsuper*a1; potStructArray.a2=Nsuper*a2;
potStructArray(1).zmin=Z(1);
potStructArray(1).zmax=Z(end);
potStructArray(1).zPoints=length(Z);
potStructArray(1).fileprefix=fileprefix;
potStructArray(1).Nxy=Nxy;
potStructArray(1).Nz=Nz;
potStructArray(1).Nsuper=Nsuper;



%% data for python hex plotter
WritePythonInfo2(fileprefix,a1,a2,b1,b2,Nsuper);

%% Get min and max bounds of the potentials
DFTmin = min(VDFTsuper,[],"all");
DFTmax = max(VDFTsuper,[],"all");
AnalyticMin = min(Vsuper,[],"all");
AnalyticMax = max(Vsuper,[],"all");

%% Now change everything to be Min's DFT
copyDFT = false;
if copyDFT
  Nsuper = DFTsuper;
  Nxy = 12;
  Nz = 19;
  Xsuper = XDFTsuper;
  Ysuper = YDFTsuper;
  Z = ZDFT;
  Vsuper = VDFTsuper;
  a1 = x1;a2=x2;
  b1 = y1; b2 = y2;
end

confStruct=Multiscat.createConfigStruct(potStructArray);
Multiscat.prepareConfigFile(confStruct);
Multiscat.prepareFourierLabels(Vsuper,fileprefix);
Multiscat.PreparePotentialFiles(potStructArray);
disp("Finished :D")

if(plotPot)
  Vplotted = Vsuper;
  
    figure
    equipotential_plot('V', Vplotted, 'V0', 0, 'z', Z, 'X', Xsuper, 'Y', Ysuper)
    shading interp
    hold on
    view([15 45])
    %equipotential_plot('V',VDFTsuper,'V0', Vsoup, 'z',ZDFT,'X',XDFTsuper,'Y',YDFTsuper)
    shading interp
    xlim([0 const.a]);
    ylim([0 const.a]);
    %xlim([-3.5 2]*Nsuper);
    %ylim([-0.5 3]*Nsuper);
    daspect([1 1 1])
    hold off
savestr = "Figures/" + fileprefix + ".jpg";
saveas(gcf,savestr,'jpg')
end
%===
%% Function definitions

function [b1,b2,b3] = Reciprocal(a1,a2,a3)
    factor = 2*pi/dot(cross(a1,a2),a3);
    b1 = factor*cross(a2,a3);
    b2 = factor*cross(a3,a1);
    b3 = factor*cross(a1,a2);
end

function [VmatrixElement] = Vfunc_LiF_Wolken(x,y,z)
    function [V0] = V0func(z)
        V0 = const.D * exp(2*const.alpha*(const.z0-z))...
            - 2*const.D*exp(const.alpha*(const.z0-z));
    end
    function [V1] = V1func(z)
        V1 = -2*const.beta*const.D*exp(2*const.alpha*(const.z0-z));
    end
    function [Q] = Qfunc(x,y)
        Q = cos(2*pi*x/const.a) + cos(2*pi*y/const.a);
    end
    VmatrixElement = V0func(z) + V1func(z)...
        * Qfunc(x,y);
end

function [VmatrixElement] = Vfunc_MoS2(X,Y,Z)
  function [V] = VSulph(z)
    D = 19.9886;
    a = 0.8122;
    alpha = 1.4477;
    b = 0.1958;
    beta = 0.2029;
    z0 = 3.3719;
    z1 = 1.7316;
    V = D*(exp(2*alpha*(z0-z))-2*a*exp(alpha*(z0-z))-2*b*exp(2*beta*(z1-z)));
  end

  function [V] = VHollow(z)
    D = 24.9674;
    a = 0.4641;
    alpha = 1.1029;
    b = 0.1993;
    beta = 0.6477;
    z0 = 3.1411;
    z1 = 3.8323;
    V = D*(exp(2*alpha*(z0-z))-2*a*exp(alpha*(z0-z))-2*b*exp(2*beta*(z1-z)));
  end

  function [V] = VMolyb(z)
    D = 20.1000;
    a = 0.9996;
    alpha = 1.1500;
    b = 0.0026;
    beta = 1.2439;
    z0 = 3.2200;
    z1 = 4.1864;
    V = D*(exp(2*alpha*(z0-z))-2*a*exp(alpha*(z0-z))-2*b*exp(2*beta*(z1-z)));
  end
        %+ V1func(Z) * Qfunc(X,Y)...
    VmatrixElement = VSulph(Z) ... %blue, sulphur
       * Qhexfunc(X,Y) ...
       + VHollow(Z) ... %green, hollow site
      * Qhexfunc(X,Y - (const.c/sqrt(3))) ...
      + VMolyb(Z) ...%red, molybdenum
      * Qhexfunc(X-const.c/2,Y-(const.c*1/(2*sqrt(3))));
end

function [Vout] = AddSulphurDefect(doWeRepeat,Vin,m_in,n_in,a1,a2,Nsuper,Xsuper,Ysuper,Z)
%Adds a defect at sulphur site (m,n)
  Vout = Vin;
  NxySuper = size(Vout,1);
  Nz = size(Vout,3);
  centre0 = double(m_in)*a1+double(n_in)*a2;
  if doWeRepeat
    disp("Repeating!")
    for m = -1:1
      for n = -1:1
        for k = 1:Nz
          centre = [centre0(1)+m*a1(1)*Nsuper+n*a2(1)*Nsuper
                    centre0(2)+m*a1(2)*Nsuper+n*a2(2)*Nsuper];
          %r = (x-centre(1))^2+(y-centre(2))^2;
          %r = sqrt(r)/const.c;
          %disp("r = " + num2str(r))
          %disp("vmatrixelem = " + Vout(i,j,k))
          %disp("defect val  = " + val(x,y,Z(k),centre))
          %disp("sum         = " + num2str(Vout(i,j,k) + val(x,y,Z(k),centre)));
          Vout(:,:,k) = Vout(:,:,k)+val(Xsuper,Ysuper,Z(k),centre);
          %disp("------------- ")
          %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
          %     ", Value = " + val(x,y,Z,k,centre));
        end
      %disp("m, n = " + m + ", " + n)
      %disp("Centre:")
      %disp(centre)
      %equipotential_plot('V', Vout, 'V0', 0, 'z', Z, 'X', Xsuper, 'Y', Ysuper)
      %hold on
      %plot3(centresX(m+2,n+2),centresY(m+2,n+2),5,'.');
      %hold off
      %xlim([-3.5 2]*Nsuper);
      %ylim([-0.5 3]*Nsuper);
      %daspect([1 1 1])
      %shading interp
      end
    end
  else
    disp("Not repeating!")
    for k = 1:Nz
        Vout(:,:,k) = Vout(:,:,k)+val(Xsuper,Ysuper,Z(k),centre0);
        %disp("x, y, z = " + x + ", " + y + ", " + Z(k) +...
        %    ", Value = " + val);
    end
  end

  function [v] = val(x,y,z,centre)
    
    function [V] = VSulph(z)
      D = 19.9886;
      a = 0.8122;
      alpha = 1.4477;
      b = 0.1958;
      beta = 0.2029;
      z0 = 3.3719;
      z1 = 1.7316;
      V = D*(exp(2*alpha*(z0-z))-2*a*exp(alpha*(z0-z))-2*b*exp(2*beta*(z1-z)));
    end
    
    r = (x-centre(1)).^2+(y-centre(2)).^2;
    r = sqrt(r)./const.c;
    %maxr = min(r,[],"all");
    %minr = min(r,[],"all");
    extent = 0.3;
    cutoff = 1;
    s = extent * const.c;
    v = 0;
    c = 0.0311;
    d = 32.8260;
    e = 16.3770;
    gamma	= 0.9209;
    lambda = 0.9204;
    z2 = 5.6012;
    z3 = 3.7072;
    %r = max(0.1,r);

    flatDefect = true;
    hexDefect = false;
    if(flatDefect)
      if(hexDefect)
        ikbT = 12.9;
        mu = 0.92;
      else
        %ikbT = 15.9;
        %mu = 0.49;

        ikbT = 4;
        mu = 0.5;
      end
      VmatrixElement = Vfunc_MoS2(x,y,z);
      %VmatrixElement = Vfunc_LiF_Wolken(x,y,z);
    else
      error("This is impossible to fit, don't use this")
      if(hexDefect)
        ikbT = 17.9;
        mu = 1.6;
      else
        ikbT = 15.9;
        mu = 0.6;
      end
      VmatrixElement = VSulph(z) * Qhexfunc(x,y);
    end
    
    %args = (y-centre(2))./(x-centre(1));
    %angle = arrayfun(@(arg) atan(arg),args);
    %angle(isnan(angle))=0;
    %cutoffR = 1/sqrt(3)*cos(pi/6)./(cos(angle-(2*pi*floor((6*angle+pi)/(2*pi)))/6));
    %cutoffR = mu .* cutoffR;

    factor = (1./( 1+exp((r-mu)*ikbT) ));
    factor = (factor.*( 1+exp((-mu)*ikbT) ));
    %factor = factor./( 1+exp((r-hardCut)*10000));
    %maxf = max(factor,[],"all");
    %minf = min(factor,[],"all");
    %if(maxf > 1. )
    %  error("Max more than 1!")
    %elseif(minf < 0.)
    %    error("Min less than 0!")
    %end
    if(isnan(factor))
      error("Nans found!")
    end

    v = (-VmatrixElement + d*(exp(2*gamma*(z2-z))-2*c*exp(gamma*(z2-z)) ...
      -2*e*exp(2*lambda*(z3-z)))).*factor;
    %disp("Maxv = " + max(v,[],"all"))
    %disp("Minv = " + min(v,[],"all"))
    %disp("1./( 1+exp((-0.45)*10) = ")
    %disp(1./( 1+exp(-0.45*10) ))
  end
end


function [Q] = Qhexfunc(X,Y)
  X_n = X ./ (const.c);
  Y_n = Y ./ (const.c*sqrt(3));
  Q = ((cos(2*pi*(X_n-Y_n))+cos(4*pi*Y_n)+cos(2*pi*(X_n+Y_n))) + 3/2)/(4.5);
  %Q = 1;
end
