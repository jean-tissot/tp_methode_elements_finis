close all ; clc;
%TP1 - partie sur les treillis
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
disp('structure etudiee : treillis du TP1');
disp('==================');

% ### VARIABLES DEPENDANT DU PROBLÈME (définition du maillage) ###
h = 1.5e3; a = 3e3; b= 2e3; % en mm
Fd=-950e3; % en N
S1=49e2; S2=25e2; % => en mm²
E=210e3; Re=350; % => en MPa
S = [S2 S1 S1];

Coord=[ 0 , 0 ;          % definition des coordonnees des noeuds X , Y
        a+b , 0 ; 
        a , -h ];   
 
Connec=[ 1 , 2 ;         % definition de la matrice de connectivite i , j
         1 , 3 ; 
         2 , 3 ];

% definition des caracteristiques mecaniques elementaires (ES fx fy)
Nprop=[2;1;1];            % pour chaque element numero de la propriete
Prop=[ E*S1 0 0 ;         % tableau des differentes valeurs de ES fx fy
       E*S2 0 0 ];

% definition des CL en deplacement
CL=[ 1 , 0 , 1 ;          % numero du noeud, type sur u et v (1 ddl impose ,0 ddl libre)
     2 , 1 , 1 ];         % + ligne 3 , 0 , 0   inutile (car aucune contrainte sur le noeud 3)

% definition des charges nodales
Charg=[ 3  0  Fd ];                %  numero du noeud , Fx , Fy

% definition du modele EF : type des elements
Typel = ['barre_ke';'barre_ke';'barre_ke']

% nombre de degrés de liberté par noeud
nddln=2;


% ### TRAITEMENT (indépendant du problème) ###

% Traitement du nombre de noeud et de la connectivité
[nnod,ndim]=size(Coord);    % nombre de noeuds / nombre de dimensions
nddlt=nddln*nnod;           % nombre de degrés de liberté total
[nelt,nnode]=size(Connec);  % nombre d'éléments / nombre de noeuds par élément

% Traitement des CL
Ncl=zeros(1,nddlt); ncld=0;  %Ncl: 1 si le de déplacement champ est imposé
Vcl=zeros(1,nddlt);          % Valeurs des deplacements imposes
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 
           Ncl(1,(CL(i,1)-1)*nddln+j)=1;
           ncld++; 
       end
   end
end

% Traitement des charges nodales
F=zeros(nddlt,1);	    %----- vecteur sollicitation
for iclf=1:size(Charg,1)
	noeud=Charg(iclf,1);
	for i=1:nddln
    F((noeud-1)*nddln+i)=F((noeud-1)*nddln+i) + Charg(iclf,i+1);
  end
end
[Fx,Fy,Fz] = feval('resultante',F);      %----- resultante des charges nodales

% Trace du maillage pour validation des donnees 
plotstr                      
U = zeros(nddlt,1);
R = zeros(nddlt,1);
[U(:,1), R(:,1)] = statique;   %----- resolution du probleme - methode du terme unite sur la diagonale    
plotdef(U)

nouv_coord = Coord + transpose(reshape(U, size(Coord,2), size(Coord,1))) %calcul des nouvelles positions des noeuds                      

%-----	post-traitement
format = [strjoin({' %8.3e',' %8.3e',' %8.3e'}(1:nddln), '  '), ' \n']; %format d'impression des vecteurs

disp(' ');
disp('------- deplacements nodaux sur (x,y,z) ----------');
fprintf(format,U)                        

disp(' ');
disp('------- Efforts aux appuis  ----------');
fprintf(format,R(:,1));

[Rx,Ry,Rz] = feval('resultante',R);     %----- resultantes et reactions
disp(' ');
fprintf('La resultante des charges nodales    en (x,y,z) est : %8.3e   %8.3e   %8.3e \n', Fx, Fy, Fz);                    
fprintf('La resultante des charges reparties  en (x,y,z) est : %8.3e   %8.3e   %8.3e \n', -Rx-Fx,-Ry-Fy, -Rz-Fz);
fprintf('La resultante des efforts aux appuis en (x,y,z) est : %8.3e   %8.3e   %8.3e \n', Rx, Ry, Rz);

disp(' ');
disp('------- Contraintes sur les elements ----------');


max_abs_contrainte = 0;
max_abs_rapport_Ne_Fc = 0;

for iel = 1:nelt          %----- boucle sur les elements
  loce=[]; 
  for i = 1:nnode
    loce=[loce, (Connec(iel,i)-1)*nddln + [1:nddln] ];
  end
  Ue = U(loce);
  Ne = feval('barre_stress',iel,Ue);

  fprintf('Dans l''element %3i l''effort normal est %8.3e\n',iel,Ne)

  num_noeud_1=Connec(iel,1); % numéro des noeuds de l'élément
  num_noeud_2=Connec(iel,2);
  % longueur de l'élément (après déformation):
  lc = norm([nouv_coord(num_noeud_1,1) nouv_coord(num_noeud_1,2)]) - norm([nouv_coord(num_noeud_2,1) nouv_coord(num_noeud_2,2)]);
  ES = Prop(Nprop(iel),1);
  % charge critique d'Euler (pi².E.I/lc² avec I=pi.r⁴/4)
  Fc = pi*ES**2/(4*E*lc**2);  % pi*I = pi².r⁴/4 = S²/4 = (E.S)²/4E
  % on garde en mémore la contrainte maximale et le rapport maximal Ne/Fc
  max_abs_rapport_Ne_Fc = max(max_abs_rapport_Ne_Fc, abs(Ne / Fc));
  max_abs_contrainte = max(max_abs_contrainte, abs(Ne * E / ES));
end

% Le calcul étant linéaire on peut en déduire la charge limite élastique et la charge limite de flambement
Fmax_elastique = Fd *  Re / max_abs_contrainte;
Fmax_euler = Fd / max_abs_rapport_Ne_Fc

fprintf('force maximale admissible quant à l''élasticité: %0.4e\n', Fmax_elastique)
fprintf('force maximale admissible quant au flambement: %0.4e\n', Fmax_euler)