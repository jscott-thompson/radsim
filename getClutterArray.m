function Clutter=getClutterArray(rt,clut_dr,clut_ring)

r_start=rt-0.5*clut_ring;
r_end=rt+0.5*clut_ring;
nr=floor(clut_ring/clut_dr);

clut_dphi=asin(clut_dr/r_end);  %   angle increment defining the clutter ring
nphi=floor(2*pi/clut_dphi);

phi=linspace(-pi/2,pi/2,nphi);
r=linspace(r_start,r_end,nr);
[Phi,R]=meshgrid(phi,r);

phig=reshape(Phi,numel(Phi),1);
rg=reshape(R,numel(R),1);

x=rg.*cos(phig);    y=rg.*sin(phig);    z=zeros(numel(rg),1);
Clutter=[x,y,z];

% figure; plot3(Clutter(:,1),Clutter(:,2),Clutter(:,3),'o');