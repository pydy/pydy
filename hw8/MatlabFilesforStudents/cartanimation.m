function []=cartanimation(qmatrix,par);
%function []=cartanimation(qmatrix,par);
%
%This function animates a cart given generalized coordinates (GCs) as
%functions of time, as specified in the qmatrix.
%These generalized coordinates of the cart are:
%2 for translation: sX, sY (coordinates of the center of mass)
%1 for rotation: theta
%The three columns of the qmatrix must contain time series of these GCs.
%parameters contain geometric descriptions of the cart.
%
%You can test this function with some fake data:
% Tstest=.01;
% t=(0:Tstest:20)';
% fakedata=[3*cos(2*pi*1*t),3*sin(2*pi*1*t),2*pi*1*t+pi/2];
% cartanimation(fakedata,par)
%
%Author: H. Vallery, October 2014

Ts_anim=par.Ts_anim;%sampling time of animation (used to pause)
numsamp=size(qmatrix,1);%number of samples

%----------------------------
%initialize the figure:
%----------------------------

%initial values:
sX=qmatrix(1,1);
sY=qmatrix(1,2);
theta=qmatrix(1,3);

rs=[sX;sY];%position of point S (center of mass) in N coordinates


%geometric data for aninmation of the cart:
length_cart=par.length_cart;%[m], length of the cart
width_cart=par.width_cart;%[m] width of the cart

%points describing the cart's geometry in body coordinates (B frame):
cart_body_B=[-length_cart/2 -length_cart/2 length_cart/2 -length_cart/2;%x coordinates of body endpoints in body coordinates
    -width_cart/2 width_cart/2 0 -width_cart/2];%y coordinates
numpbody=size(cart_body_B,2);%number of points


%rotation matrix construction:

R_theta=[cos(theta) sin(theta);
    -sin(theta) cos(theta)];%rotation about z axis



%geometric data expressed in the N frame:
cart_body_N= R_theta'*cart_body_B+repmat(rs,1,numpbody);
track_N= cart_body_N(:,1:2);

%visualize:
figure(88)
clf
handle_cart_body=plot(cart_body_N(1,:),cart_body_N(2,:),'b','linewidth',4);
hold all
handle_cart_wheels=plot(cart_body_N(1,:),cart_body_N(2,:),'r*','linewidth',4);

plot(0,0,'bx')


xlabel('x')
ylabel('y')



axis equal
maxdim=5;
xlim([-maxdim,maxdim])
ylim([-maxdim,maxdim])

%----------------------------
%animate by iterating through all following samples:
%----------------------------

for indexsamples=1:numsamp
    sX=qmatrix(indexsamples,1);
    sY=qmatrix(indexsamples,2);
    theta=qmatrix(indexsamples,3);
    
    rs=[sX;sY];%position of point S in N coordinates
    
    R_theta=[cos(theta) sin(theta);
        -sin(theta) cos(theta)];%rotation about z axis
   
    
    %geometric data expressed in the N frame:
    cart_body_N= R_theta'*cart_body_B+repmat(rs,1,numpbody);
    track_N= cart_body_N(:,1:2);
 
    
    %update the data in the plot (faster than new plot each time):
    set(handle_cart_body,'XData',cart_body_N(1,:),'YData',cart_body_N(2,:));
    set(handle_cart_wheels,'XData',cart_body_N(1,:),'YData',cart_body_N(2,:));
    plot(track_N(1,:),track_N(2,:),'k.');

    pause(Ts_anim)%generate soft realtime
    
    
    
%    cd film_cart %DIESES VERZEICHNIS MUSS SCHON EXISTIEREN
%   % set(gcf,'position',[10 50 640 480]);
%   set(gcf,'position',[500 500 280 280]);
%    set(gcf,'color','w');
%    eval(['print -f' num2str(gcf) ' -dbitmap cartanimation_' num2str(indexsamples) '.bmp'])
%    cd ..
% %also need to change Ts in simulation to .04, to get 25
% %frames per second for a video

end

