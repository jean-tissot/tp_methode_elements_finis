close all ; clc;
%  jeu de donnees du treillis traite en exemple dans le cours
%  avec a=100 , ES = 100*sqrt(2) , F = 40.
%  avec la methode du terme unite sur la diagonale
%
% Fonctions utilisees
%   statique        : calcul de la reponse statique [U,R]
%   plotstr         : trace du maillage avec ne noeuds et elements
%   plotdef         : trace de la deformee
%   resultante      : calcul de la resultante en x, y, z d'un vecteur nodal
%   barre_stress    : calcul de la contrainte dans un element barre

global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F
disp(' ');
disp('structure etudiee : treillis traite en exemple dans le cours');
disp('==================');
% definition du maillage
h = 1.5e3; a = 3e3; b= 2e3; S1=49e2; S2=25e2; E=210e3 %unités en mm et N (=> mm² et MPa)
S = [S2 S1 S1]
Re=350 %MPa
Fd=-950e3
Coord=[ 0 , 0 ; ...         % definition des coordonnees des noeuds X , Y
        a+b , 0 ; ...
        a , -h ];
[nnod,ndim]=size(Coord);  %nombre de noeuds / nombre de dimensions (3,2)
nddln=2;  nddlt=nddln*nnod; %nombre de degrés de liberté par noeud / nbre total    
 
Connec=[ 1 , 2 ; ...        % definition de la matrice de connectivite i , j
         1 , 3 ; ...
         2 , 3 ];
[nelt,nnode]=size(Connec); %nombre d'éléments, nombre de noeuds par élément (3,2)

%definition du modele EF : type des elements
Typel = 'barre_ke';               % definition du type des elements
for i=1:nelt
    Typel = char('barre_ke',Typel);
end
% definition des caracteristiques mecaniques elementaires (ES fx fy)
Nprop=[2;1;1];              % pour chaque element numero de la propriete
Prop=[ E*S1 0 0 ;     % tableau des differentes valeurs de ES fx fy
      E*S2 0 0 ;
      ];
      
% definition des CL en deplacement
CL=[ 1 , 0 , 1 ; ... % numero du noeud, type sur u et v (1 ddl impose ,0 ddl libre)
     2 , 1 , 1 ];  % + ligne 3 , 0 , 0   inutile (car aucune contrainte sur le noeud 3)
Ncl=zeros(1,nddlt); ncld=0;  %Ncl: 1 si le de déplacement champ est imposé
Vcl=zeros(1,nddlt);       % Valeurs des deplacements imposes
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 
           Ncl(1,(CL(i,1)-1)*nddln+j)=1;
           ncld=ncld+1; 
       end
   end
end

% definition des charges nodales
Charg=[ 3  0  Fd ];                %  numero du noeud , Fx , Fy
      
F=zeros(nddlt,1);	    %----- vecteur sollicitation
for iclf=1:size(Charg,1)           
	noeud=Charg(iclf,1);  
	for i=1:nddln
    F((noeud-1)*nddln+i)=F((noeud-1)*nddln+i) + Charg(iclf,i+1);
    end
    end
[Fx,Fy,Fz] = feval('resultante',F);      %----- resultante des charges nodales

% trace du maillage pour validation des donnees 
close all
plotstr                      
U = zeros(nddlt,1);
R = zeros(nddlt,1);
[U(:,1),R(:,1)] = statique;   % ----- resolution du probleme
%                               methode du terme unite sur la diagonale                              
%----- format d'impression des vecteurs
form =' %8.3e   %8.3e   %8.3e  '; format = [form(1:8*nddln),' \n']; 
disp(' ');disp('------- deplacements nodaux sur (x,y,z) ----------');
fprintf(format,U)
plotdef(U)                                             
%-----	post-traitement
disp(' ');disp('------- Efforts aux appuis  ----------');
fprintf(format,R(:,1));
[Rx,Ry,Rz] = feval('resultante',R);     %----- resultantes et reactions

disp(' ');
fprintf('La resultante des charges nodales    en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',Fx,Fy,Fz);                    
fprintf('La resultante des charges reparties  en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',-Rx-Fx,-Ry-Fy,-Rz-Fz);
fprintf('La resultante des efforts aux appuis en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',Rx,Ry,Rz);

disp(' ');disp('------- Contraintes sur les elements ----------');
max_abs_contrainte=0;
liste_Ne=[0 0 0]
for iel=1:nelt          %----- boucle sur les elements
  loce=[]; for i=1:nnode loce=[loce,(Connec(iel,i)-1)*nddln+[1:nddln]];end
  Ue=U(loce);
  Ne=feval('barre_stress',iel,Ue);
  contrainte_abs=abs(Ne/S(iel));
  max_abs_contrainte=max(max_abs_contrainte, contrainte_abs)
  liste_Ne(iel)=Ne
  fprintf('Dans l''element %3i l''effort normal est %8.3e\n',iel,Ne)
end

Fmax = Fd *  Re / max_abs_contrainte; %calcul de la limite

fprintf('force maximale admissible: %0.4e\n', Fmax)

%calcul des nouvelles positions (et des longueurs)
nouv_coord = Coord + transpose(reshape(U,2,3))
longueurs=[0 0 0];
for i=1:3
  num_noeud_1=Connec(i,1);
  num_noeud_2=Connec(i,2);
  longueurs(i)=norm([nouv_coord(num_noeud_1,1) nouv_coord(num_noeud_1,2)]) - norm([nouv_coord(num_noeud_2,1) nouv_coord(num_noeud_2,2)]);
end
disp("longueur des 3 barres après déformation:")
disp(longueurs)

liste_Fc=[0 0 0];
for i=1:3
  lc=longueurs(i);
  liste_Fc(i)=pi*E*S(i)**2/(4*lc**2);  % pi*I = pi²*r**4/4 = S²/4
end

disp("Charge critique pour les 3 barres:")
disp(liste_Fc)

rapport=min(abs(liste_Fc./liste_Ne));
Fc=Fd*rapport;

disp("Charge critique de flambement:")
disp(Fc);