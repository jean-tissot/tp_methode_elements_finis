close all ; clc;
%  Script pour la poutre de l'exercice de cours exo15
%         poutre sur 2 appuis, modelisee par nelt elements
% 
%  Fonctions utilisees
%   statiqueUR      : calcul de la reponse statique [U,R]
%   plotstr         : trace du maillage avec numero noeuds et elements
%   plotdef         : trace de la deformee
%   resultante      : calcul de la resultante en x, y, z d'un vecteur nodal
%   poutre_stress   : calcul de la contrainte dans un element poutre
%   poutre_compar   : comparaison avec la solution analytique
%   
global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F

disp(' ');
disp('structure etudiee : poutre de l''exercice de cours exo15');
disp('==================');

% ### VARIABLES DEPENDANT DU PROBLÈME (définition du maillage) ###
nelt=20 % nombre d'élélément de découpe de la poutre -> doit être un nombre pair

l=0.5       % demi-longueur de la poutre (en m)
f=10e3     % force linéique en N/m
EI=1

Coord=[]; 
for j=0:nelt
  Coord=[Coord; j*2*l/nelt];
end

Connec=[];
for j=1:nelt
  Connec=[Connec;[j j+1]];    % chaque element j est relié à l'élément j-1 et à l'élément j+1
end

% definition des caracteristiques mecaniques elementaires (EI f)  (en 1D)
Nprop = ones(1,nelt);           % pour chaque element numero de la propriete
Prop=[EI f]                   % -- Propriété 1 pour la deuxième moitié de la poutre

for i = 1:nelt/2
  x1=Coord(i, 1);             % abscisse gauche de l'élément i
  x2=Coord(i+1, 1);           % abscisse droite de l'élément i
  f_i = (f*x1/l + f*x2/l)/2;  % moyenne de la force linéique sur l'élément i

  Prop = [Prop; EI f_i];      % -- Propriété i+1 pour l'élément i de la première moitié de la poutre
  Nprop(i) = i+1;
end

% definition des CL en deplacement
CL=[ 1 ,        1 , 1 ;       % numero du noeud, type sur v et theta (1 ddl impose ,0 ddl libre)
     1+nelt/2 , 1 , 0 ;
     1+nelt ,   1 , 0 ];

%Charges nodales: aucune
charg=[] 


% definition du modele EF : type des elements
Typel = 'poutre_ke';    
for i = 1:nelt-1
  Typel = char('poutre_ke',Typel);
end

% nombre de degrés de liberté par noeud
nddln=2;


% ### TRAITEMENT (indépendant du problème) ###

% Traitement du nombre de noeud et de la connectivité
[nnod,ndim]=size(Coord);
nddlt=nddln*nnod;
nnode = 2; % nelt, nnode = size(Connec)  - nnode=2 car chaque poutre a 2 noeuds

% Traitement des CL
Ncl=zeros(1,nddlt);ncld=0;
Vcl=zeros(1,nddlt);         % Valeurs imposees nulles
for i=1:size(CL,1)
    for j=1:nddln 
      if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1; ncld++; end
    end
end


% definition des charges nodales
F=zeros(nddlt,1);	   
[Fx,Fy,Fz] = feval('resultante',F); 

% trace du maillage pour validation des donnees
plotstr
% ----- resolution du probleme
U = zeros(nddlt,1);
R = zeros(nddlt,1);
[U(:,1),R(:,1)] = statiqueUR;   
plotdef(U)

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
fprintf('La resultante des charges nodales    en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',Fx,Fy,Fz);                    
fprintf('La resultante des charges reparties  en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',-Rx-Fx,-Ry-Fy,-Rz-Fz);
fprintf('La resultante des efforts aux appuis en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',Rx,Ry,Rz);

disp(' ');
disp('------- Contraintes sur les elements ----------');

for iel=1:nelt          %----- boucle sur les elements
  loce=[];
  for i = 1:nnode
    loce=[loce, (Connec(iel,i)-1)*nddln + [1:nddln] ];
  end
  Ue=U(loce);
  feval('poutre_stress',iel,Ue);
end

feval('poutre_compar',U); %----- comparaison avec la solution analytique                                     
return