clc
clear all
close all
rng('default');
rng(1);

% Spatial synthesis of SS dyads

%Generate unit dual quaternions
npose=7;
AxisAng=rand(npose,4);
AxisAng(:,1:3)=AxisAng(:,1:3)./vecnorm(AxisAng(:,1:3),2,2);
q1=(AxisAng(:,1).*sin(AxisAng(:,4)/2))';
q2=(AxisAng(:,2).*sin(AxisAng(:,4)/2))';
q3=(AxisAng(:,3).*sin(AxisAng(:,4)/2))';
q4=(cos(AxisAng(:,4)/2))';

Disp=rand(npose,3);
for i=1:npose
    d1=Disp(i,1); d2=Disp(i,2); d3=Disp(i,3);
    g=.5*[0,-d3,d2,d1;d3,0,-d1,d2;-d2,d1,0,d3;-d1,-d2,-d3,0]*[q1(i);q2(i);q3(i);q4(i)];
    g1(i)=g(1);
    g2(i)=g(2);
    g3(i)=g(3);
    g4(i)=g(4);
end

% Creating the K matrix
K0=ones(1,npose);
K1=2*(q4.^2+q1.^2-q2.^2-q3.^2);
K2=4*(q1.*q2+q4.*q3);
K3=4*(q1.*q3-q4.*q2);
K4=4*(g4.*q1+g3.*q2-g2.*q3-g1.*q4);
K5=4*(q1.*q2-q4.*q3);
K6=2*(q4.^2-q1.^2+q2.^2-q3.^2);
K7=4*(q2.*q3+q4.*q1);
K8=4*(-g3.*q1+g4.*q2+g1.*q3-g2.*q4);
K9=4*(q1.*q3+q4.*q2);
K10=4*(q2.*q3-q4.*q1);
K11=2*(q4.^2-q1.^2-q2.^2+q3.^2);
K12=4*(g2.*q1-g1.*q2+g4.*q3-g3.*q4);
K13=-4*(g4.*q1-g1.*q4+g2.*q3-g3.*q2);
K14=-4*(g4.*q2-g2.*q4+g3.*q1-g1.*q3);
K15=-4*(g4.*q3-g3.*q4+g1.*q2-g2.*q1);
K16=-4*(g1.^2+g2.^2+g3.^2+g4.^2);

K=[K0',K1',K2',K3',K4',K5',K6',K7',K8',K9',K10',...
    K11',K12',K13',K14',K15',K16'];

[U,S,V] = svd(K);
nullSp=V(:,8:end);

W1=nullSp(1:4,:);
W2=nullSp(5:8,:);
W3=nullSp(9:12,:);
W4=nullSp(13:16,:);

syms lamda1 lamda2 lamda3
mat1=W2-lamda1*W1;
mat2=W3-lamda2*W1;
mat3=W4-lamda3*W1;
mat4=cat(1,mat1,mat2,mat3);
%vpa(mat4,2)
vpa(det(mat4(1:10,:)),2)
vpa(det(mat4(2:11,:)),2)
vpa(det(mat4(3:12,:)),2)

eq=[det(mat4(1:10,:))==0,det(mat4(2:11,:))==0,det(mat4(3:12,:))==0];
vars=[lamda1,lamda2,lamda3];
%sol=solve(eq,vars);