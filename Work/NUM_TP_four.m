close all; clear all; clc;
%  four :  Exemple de creation des donnees pour le cas test du four
%  pour répondre aux questions du texte de TP Therm 
%  vous pouvez partir de ce script et le modifier
%  fonctions utilisees
%        statiqueU       : resolution du Pb thermique par EF 
%        plot_therm      : carte des temperatures sur le domaine
% H. Oudin
global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F

% donnees du probleme
L=2 ; h=2 ; T1=0 ; T2=0 ; r=1 ; k=1 ;
disp(' ');
disp('Programme thermique : cas test four simple');
disp('====================');
% definition du maillage
nex = input('donner le nombre d''element en x (longueur) ? [4]: ');
if isempty(nex) nex=4; end  
ney = input('donner le nombre d''element en y (hauteur) ? [4]: ');
if isempty(ney) ney=4; end  
Coord=[];         
for i=0:nex
    for j=0:ney Coord=[Coord;[i*L/nex  j*h/ney]];end 
    end
Connec=[]; pas=ney+1;   
for i=1:nex
   for j=1:ney i1=(i-1)*pas+j; Connec=[Connec;[i1  i1+pas i1+pas+1 i1+1]]; end 
   end
[nnod,ndim]=size(Coord);
[nelt,nnode]=size(Connec);
nddln=1;  nddlt=nddln*nnod;               
% definition du type des elements
Typel = [];                     
for i=1:nelt
    Typel = [Typel; 'Q4_th'] ;         
end
% definition des caracteristiques thermiques elementaires (k r)
Nprop=ones(nex*ney);
Prop=[ k r ];           % tableau des differentes valeurs  k  r
% definition des CL en temperature       
CL=[]; Vcl=zeros(1,nddlt);     
for i=1:nex+1
   n=i*(ney+1);
   CL=[CL;[n-ney , 1]];
   CL=[CL;[n , 1]];
   Vcl(n-ney)=T1;
   Vcl(n)=T2;
end
Ncl=zeros(1,nddlt);ncld=0;
for i=1:size(CL,1)
   for j=1:nddln 
   if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1;ncld=ncld+1; end
   end
end
% definition des flux nodaux                      
Charg=[ ];              
F=zeros(nddlt,1);	    % vecteur sollicitation

U = zeros(nddlt,1);
U = statiqueU;             
taille = get(0,'ScreenSize'); 
figure('Name','Temperatures dans le four',...
        'Position',[taille(3)/2.01 taille(4)/2.6 taille(3)/2 taille(4)/2])
plot_therm(U);
return