clc 
clear all
close all

% Calculate unit dual quaternion 
% Qr= Rotation quaternion
% Qd= Translation quat

% Planar displacement
% x,y,ang
d1=10; d2=20; d3=0;
ang=pi/3;
axis=[0,0,1];
Qr=[cos(ang/2),sin(ang/2)*axis]
Trans=[0 -d3 d2 d1
    d3 0 -d1 d2
    -d2 d1 0 d3
    -d1 -d2 -d3 0];
Qd= (.5*Trans*[Qr(2:4),Qr(1)]')'

%Verifying with planar quaternion
z1=.5*(d1*cos(ang/2)+d2*sin(ang/2));
z2=.5*(d1*-sin(ang/2)+d2*cos(ang/2));
z3=sin(ang/2);
z4=cos(ang/2);
pl=[z1,z2,z3,z4];

% Spherical displacement
axis=[1,1,1];
ang=pi/3;
Qr=[cos(ang/2),sin(ang/2)*axis];