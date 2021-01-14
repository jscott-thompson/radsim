function [Array] = getArrayElements(Ny,Nz,dyz)

y=linspace(1,Ny,Ny)*dyz;
z=linspace(1,Nz,Nz)*dyz;

[Y,Z]=meshgrid(y,z);

yg=reshape(Y,Nz*Ny,1);  zg=reshape(Z,Nz*Ny,1);  xg=zeros(Nz*Ny,1);
Array=[xg,yg,zg];



