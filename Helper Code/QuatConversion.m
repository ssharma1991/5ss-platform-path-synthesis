x=[-0.32,0.4322,-0.8431];
y=[0.3764,-0.7587,-0.5318];
z=cross(x,y)


%-0.869504
%-0.487519
%0.0801039

rotm=[x',y',z'];

ax=rotm2axang(rotm)
q=axang2quat(ax)