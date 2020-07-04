% Simulate a 5ss-mechanism
clc
clear all
close all

INDEX=4480;

Pts=[[-7.20078633 -5.2930031   7.87116744]
 [ 5.71803324 -0.71371165 -8.25749538]
 [-8.74711427 -4.59885719  5.48763714]
 [-9.99593062 -8.72040051  6.08658925]
 [-3.87891519  6.61366631  8.91117982]
 [-6.60806082 -9.9868523   3.93482383]
 [ 8.19190297 -9.04651467  6.06662798]
 [ 7.83840293 -5.73238752 -0.61694473]
 [ 4.14722147 -0.48122018 -2.04096393]
 [-0.97390133 -9.04023547  1.92449292]
 [ 2.19533097 -7.1743553   5.35866566]];

[cplrpath,cplrorient]=simulate5SS(Pts);
cla
drawMech(Pts,cplrpath,cplrorient)
length(cplrpath)

function [cplr, cplrO] = simulate5SS(Pts)
[cplrPath_P,cplrOrient_P]=simulate5SS_linearActuation(Pts, .01);
[cplrPath_N,cplrOrient_N]=simulate5SS_linearActuation(Pts, -.01);
cplr=cat(1,flip(cplrPath_N),cplrPath_P(2:end,:));
cplrO=cat(1,flip(cplrOrient_N),cplrOrient_P(2:end,:));
end
function [cplr_Path, cplr_Orient] = simulate5SS_linearActuation(Pts, iter)
DistConsts=[1,6; 2,7; 3,8; 4,9; 5,10; 
    11,6; 11,7; 11,8; 11,9; 11,10; 
    6,7; 6,8; 6,9; 6,10; 
    7,8; 7,9; 7,10; 
    1,7];
P0=cat(2,Pts(6:end,:),ones(6,1))';

L_init=LenConst(Pts,DistConsts);
L_target=L_init;
disp=0;
cplr_Path=[];
cplr_Orient=[];
while(true)
    L_target(end)=L_init(end)+disp;
    Pts=NewtonRhapson(L_target,Pts,DistConsts);
    
    F=Cost(LenConst(Pts,DistConsts),L_target);
    if norm(F)>10^-8
        break
    end
    disp=disp+iter;
    cplr_Path=[cplr_Path;Pts(11,:)];
    
    % Calculate displacement matrix
    Pf=cat(2,Pts(6:end,:),ones(6,1))';
    Disp=Pf*pinv(P0);
    Q=calcRot2Quat(Disp(1:3,1:3));
    cplr_Orient=[cplr_Orient;Q];
    % Check Result : Pf-Disp*P0=0
    
    %cla
    %drawMech(Pts,cplr_Path,cplr_Orient);
    %drawnow
end

end


% SIMULATION Functions
function [Pts]= NewtonRhapson(L_target,Pts,DistConsts)
F=Cost(LenConst(Pts,DistConsts),L_target);
for i=1:10
    F_dash=CostGradient(Pts,DistConsts);
    dX=-F_dash\F;
    dX=reshape(dX,[3,6])';
    Pts(6:11,:)=Pts(6:11,:)+dX;
    F=Cost(LenConst(Pts,DistConsts),L_target);
    if (norm(F)<10^-8)
        break
    end
end
end
function [Cst]= Cost(CurrentL,TargetL)
Cst=(CurrentL-TargetL);
end
function [grad]= CostGradient(Pts,DistConsts)
grad=zeros(18,18);
for i=1:18
    % Permute over each constraint
    Pt1_index=DistConsts(i,1);
    Pt2_index=DistConsts(i,2);
    Pt1=Pts(Pt1_index,:);
    Pt2=Pts(Pt2_index,:);
    
    %Partial derivative wrt x,y,z of first coordinate
    if Pt1_index>5
        dx=(Pt1(1)-Pt2(1))/norm(Pt1-Pt2);
        dy=(Pt1(2)-Pt2(2))/norm(Pt1-Pt2);
        dz=(Pt1(3)-Pt2(3))/norm(Pt1-Pt2);
        
        grad(i,3*(Pt1_index-5)-2)=dx;
        grad(i,3*(Pt1_index-5)-1)=dy;
        grad(i,3*(Pt1_index-5))=dz;
    end
    
    %Partial derivative wrt x,y,z of second coordinate
    dx=-(Pt1(1)-Pt2(1))/norm(Pt1-Pt2);
    dy=-(Pt1(2)-Pt2(2))/norm(Pt1-Pt2);
    dz=-(Pt1(3)-Pt2(3))/norm(Pt1-Pt2);
        
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
Lengths=Lengths';
end

% DRAWING function
function []= drawMech(Pts,cplr_Path,cplr_Orient)
dyads=[Pts(1:5,:),Pts(6:10,:)];
cplr=Pts(11,:);

plot3(0,0,0,'*')
hold on

%Print SS Dyads
for i=1:5
    x=[dyads(i,1),dyads(i,4)];
    y=[dyads(i,2),dyads(i,5)];
    z=[dyads(i,3),dyads(i,6)];
    plot3(x,y,z,'LineWidth',2,'Color',[.5 .5 .5]);
end

%Print Coupler
for i=1:5
    x=[dyads(i,4),cplr(1)];
    y=[dyads(i,5),cplr(2)];
    z=[dyads(i,6),cplr(3)];
    plot3(x,y,z,'LineWidth',2,'Color','#0072BD')
end

%Print Actuation link
x=[dyads(1,1),dyads(2,4)];
y=[dyads(1,2),dyads(2,5)];
z=[dyads(1,3),dyads(2,6)];
plot3(x,y,z,':b','LineWidth',2);

%Print Fixed pivots
for i=1:5
    scatter3(dyads(i,1),dyads(i,2),dyads(i,3),8,'k^','LineWidth',4)
end

%Print Moving pivots
for i=1:5
    scatter3(dyads(i,4),dyads(i,5),dyads(i,6),2,'o','LineWidth',8,'MarkerEdgeColor','#0072BD')
end

%Print Coupler point
scatter3(cplr(1),cplr(2),cplr(3),20,'ko','LineWidth',1)
axis ([-11 11 -11 11 -11 11])
pbaspect([1 1 1])

%Print Coupler Path
plot3(cplr_Path(:,1),cplr_Path(:,2),cplr_Path(:,3),'LineWidth',2,'Color',[.7 .7 .7]);

%Print Coupler Orientation
[n,~]=size(cplr_Orient);
for i=1:round(n/25):n
    Ci=cplr_Path(i,:);
    Qi=cplr_Orient(i,:);
    drawCoordSys(Qi,Ci)
end

end
function []= drawCoordSys(Q,Origin)
rot=quat2rotm([Q(4),Q(1),Q(2),Q(3)])';
O=Origin;
%x axis
quiver3(O(1),O(2),O(3),rot(1,1),rot(2,1),rot(3,1),1,'LineWidth',1,'Color','r','MaxHeadSize',1)
%y axis
quiver3(O(1),O(2),O(3),rot(1,2),rot(2,2),rot(3,2),1,'LineWidth',1,'Color',[ 0.4660,0.6740,0.1880],'MaxHeadSize',1)
%z axis
quiver3(O(1),O(2),O(3),rot(1,3),rot(2,3),rot(3,3),1,'LineWidth',1,'Color','b','MaxHeadSize',1)
end


% Helper functions
function [Q]= calcRot2Quat(Rot)
Qi=rotm2quat(Rot);
Q=[Qi(2),Qi(3),Qi(4),Qi(1)];
end
