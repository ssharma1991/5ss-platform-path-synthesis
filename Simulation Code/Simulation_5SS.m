clc
clear all
close all

%Dyads matrix [Fx,Fy,Fz,Mx,My,Mz]
rng(1);
dyads=random('Uniform',-10,10,[5,6]);
cplr=random('Uniform',-10,10,[1,3]);
Pts=[dyads(1:5,1:3);dyads(1:5,4:6);cplr];
simulate5SS(Pts)

function [] = simulate5SS(Pts)
drawMech(Pts);
DistConsts=[1,6; 2,7; 3,8; 4,9; 5,10; 
    11,6; 11,7; 11,8; 11,9; 11,10; 
    6,7; 6,8; 6,9; 6,10; 
    7,8; 7,9; 7,10; 
    1,7];
L=LenConst(Pts,DistConsts);
F=Cost(L,L)
F_dash=CostGradient(Pts,DistConsts,L);
inv(F_dash)
end


% SIMULATION Functions
function [Cst]= Cost(initLen,newLen)
Cst=(newLen-initLen);
end
function [grad]= CostGradient(Pts,DistConsts,initLen)
grad=zeros(18,18);
for i=1:18
    % Permute over each constraint
    Pt1_index=DistConsts(i,1);
    Pt2_index=DistConsts(i,2);
    Pt1=Pts(Pt1_index,:);
    Pt2=Pts(Pt2_index,:);
    
    %Partial derivative wrt x,y,z of first coordinate
    if Pt1_index>5
        dx=(Pt1(1)-Pt2(1))/initLen(i);
        dy=(Pt1(2)-Pt2(2))/initLen(i);
        dz=(Pt1(3)-Pt2(3))/initLen(i);
        
        grad(i,3*(Pt1_index-5)-2)=dx;
        grad(i,3*(Pt1_index-5)-1)=dy;
        grad(i,3*(Pt1_index-5))=dz;
    end
    
    %Partial derivative wrt x,y,z of second coordinate
    dx=-(Pt1(1)-Pt2(1))/initLen(i);
    dy=-(Pt1(2)-Pt2(2))/initLen(i);
    dz=-(Pt1(3)-Pt2(3))/initLen(i);
        
    grad(i,3*(Pt2_index-5)-2)=dx;
    grad(i,3*(Pt2_index-5)-1)=dy;
    grad(i,3*(Pt2_index-5))=dz;
end
end
function [Lengths]= LenConst(Pts,DistConsts)
n_const=length(DistConsts);
for i=1:n_const
    Pt1_index=DistConsts(i,1);
    Pt2_index=DistConsts(i,2);
    Pt1=Pts(Pt1_index,:);
    Pt2=Pts(Pt2_index,:);
    Lengths(i)=norm(Pt1-Pt2);
end
end

% DRAWING function
function []= drawMech(Pts)
dyads=[Pts(1:5,:),Pts(6:10,:)];
cplr=Pts(11,:);

plot3(0,0,0,'*')
hold on

%Print SS Dyads
for i=1:5
    x=[dyads(i,1),dyads(i,4)];
    y=[dyads(i,2),dyads(i,5)];
    z=[dyads(i,3),dyads(i,6)];
    plot3(x,y,z,'b','LineWidth',2);
end

%Print Coupler
for i=1:5
    x=[dyads(i,4),cplr(1)];
    y=[dyads(i,5),cplr(2)];
    z=[dyads(i,6),cplr(3)];
    plot3(x,y,z,'g','LineWidth',2)
end

%Print Actuation link
x=[dyads(1,1),dyads(2,4)];
y=[dyads(1,2),dyads(2,5)];
z=[dyads(1,3),dyads(2,6)];
plot3(x,y,z,':b','LineWidth',2);

%Print Fixed pivots
for i=1:5
    scatter3(dyads(i,1),dyads(i,2),dyads(i,3),8,'k^','LineWidth',2)
end

%Print Moving pivots
for i=1:5
    scatter3(dyads(i,4),dyads(i,5),dyads(i,6),2,'bo','LineWidth',2)
end

%Print Coupler point
scatter3(cplr(1),cplr(2),cplr(3),20,'ko','LineWidth',1)
axis ([-10 10 -10 10 -10 10])
pbaspect([1 1 1])
end
