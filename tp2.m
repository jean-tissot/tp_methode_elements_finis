%  NUM_2
%  Poutre console avec un mailllage automatique en Q4
%  H.Oudin 

global nddln nnod nddlt nelt nnode ndim ncld
global Coord Connec Typel Nprop Prop Ncl Vcl F
clc;
L = 3000 ; h =200 ; e=10;    % donnees unites N , mm
E=210.e+3; nu=0.3;
Mf=4e3;
poid = 100;

%nex=8; ney=4;                  % definition du maillage
nex = input('donner le nombre d''element en x (longueur) ? [5]: ');
if isempty(nex) nex=5; end  
ney = input('donner le nombre d''element en y (hauteur) ? [4]: ');
if isempty(ney) ney=4; end  

Coord=[];                                               %maillage
    for i=0:nex
        for j=0:ney Coord=[Coord;[i*L/nex  j*h/ney]];end 
    end
Connec=[]; pas=ney+1;   
    for i=1:nex
        for j=1:ney 
            i1=(i-1)*pas+j;
            Connec=[Connec;[i1  i1+pas i1+pas+1 i1+1]];
        end 
    end
[nnod,ndim]=size(Coord);
[nelt,nnode]=size(Connec);
nddln=2;  nddlt=nddln*nnod; 
Typel = [];                
for i=1:nelt
    Typel = char('Q4_ep',Typel);
end

Nprop=ones(nex*ney);   % caract�ristiques m�caniques (E nu e fx fy)
Prop=[ E nu e 0 0 ];

CL=[1 1 1; nnod-ney 0 1];          % definition des CL en d�placement

Ncl=zeros(1,nddlt);ncld=0;
Vcl=zeros(1,nddlt);        
for i=1:size(CL,1)
   for j=1:nddln 
       if CL(i,1+j)==1 Ncl(1,(CL(i,1)-1)*nddln+j)=1;ncld=ncld+1; end
   end
end

Charg=[nnod , -Mf/h, 0 ; nnod-ney, Mf/h, 0];     % definition des charges nodales
% somme_bras_de_levier = 0;
% bras_de_levier_elementaire = h/ney;
% for i=1:ney/2
%     somme_bras_de_levier += 2*bras_de_levier_elementaire*i;
% end
% f_lineique = (Mf/somme_bras_de_levier)/h;
% Charg = [];
% for i=0:ney-1
%     Charg = [Charg; [nnod-ney+i, f_lineique*(ney/2 - i)*bras_de_levier_elementaire, 0]];
% end

F=zeros(nddlt,1);	 
for iclf=1:size(Charg,1)           
	noeud=Charg(iclf,1);  
	for i=1:nddln
       F((noeud-1)*nddln+i)=F((noeud-1)*nddln+i) + Charg(iclf,i+1);
    end
 end
[Fx,Fy,Fz] = feval('resultante',F);      %----- r�sultante des charges nodales

%clc;
disp(' ');
disp('Les variables globales sont initialisees');
disp('Fin de lecture des donnees');
close all
plotstr  % trac� du maillage pour validation des donn�es 
disp(' ');
disp('Elasticite plane  : Plaque en contrainte plane modelisee en Q4');
disp('==================');
U = zeros(nddlt,1);
R = zeros(nddlt,1);
[U(:,1),R(:,1)] = statiqueUR;   % ----- r�solution du probl�me

                            %----- format d'impression des vecteurs
form =' %8.3e   %8.3e   %8.3e  '; format = [form(1:8*nddln),' \n']; 
disp(' ');disp('------- deplacements nodaux sur (x,y,z) ----------');
%fprintf(format,U)
plotdef(U)
                                        %-----	post-traitement
disp(' ');disp('------- Efforts aux appuis  ----------');
%fprintf(format,R(:,1));
[Rx,Ry,Rz] = feval('resultante',R);     %----- r�sultantes et r�actions

disp(' ');
fprintf('La resultante des charges nodales    en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',Fx,Fy,Fz);                    
fprintf('La resultante des charges reparties  en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',-Rx-Fx,-Ry-Fy,-Rz-Fz);
fprintf('La resultante des efforts aux appuis en (x,y,z) est : %8.3e   %8.3e   %8.3e \n',Rx,Ry,Rz);


disp(' ');disp('------- Contraintes sur les elements ----------');
sige=[];   % tableau des contraintes sigxx
for iel=1:nelt          %----- boucle sur les �l�ments
    loce=[]; 
    for i=1:nnode 
		if Connec(iel,i) > 0 loce=[loce,(Connec(iel,i)-1)*nddln+[1:nddln]]; end
	end;  
	Ue=U(loce);
    sig = feval([deblank(Typel(iel,:)),'_stress'],iel,Ue);
    sVM = sqrt(sig(1)^2+sig(2)^2-sig(1)*sig(2)+3*sig(3)^2);   % contrainte moyenne de Von Mises sur l'�l�ment
    VM(iel)=sVM ;
%    fprintf('Valeur Moyenne des contraintes dans l''element %3i\n',iel)
%    fprintf('sigxx = %8.3e  sigyy = %8.3e   sigxy = %8.3e  sigVM = %8.3e \n',sig(1),sig(2),sig(3),sVM)
    sige=[sige;[sig(1)]];
end
figure('Name',"contraintes")
plot_sig(VM)  

% comparaison avec la solution poutre
I=e*h^3/12; EI=E*I ;        % solution poutre
v=poid*(L^3)/(3*EI); sigxx=poid*L*h/(2*I);
disp(' ');
fprintf('Solution analytique poutre\n')
fprintf('la fleche en bout est = %8.3e mm et la contrainte max sigxx = %8.3e  MPa \n',v,sigxx)
disp(' ');

%solution analytique
%Ei*v,x² = Mf.x/l
%Ei*v = Mf.x³/(6l) + a.x + b
%v(O) = 0 => b = 0
%v(l) = 0 => Mf.l²/6 + a.l = 0 => a = -Mf.l/6
%v = (Mf/(6.Ei.l)).x³ - (Mf.l/(6.Ei)).x
%v = (Mf/(6.Ei))*[x*[(1/l).x²-l]]

T=linspace(0,L,100);
for i=1:size(T,2)
    Sol(i)=Mf*(T(i)*((T(i)^2)/L-L))/(6*EI);
end
taille = get(0,'ScreenSize'); 
figure('Name','comparaison des fleches 2D-Q4 (rouge) / poutre (bleu)')  
plot(T,Sol);
hold on
Uv=U([2:2*(ney+1):2*(nex+1)*(ney+1)]);
Tv=linspace(0,L,size(Uv,1));
plot(Tv,Uv,'r');

figure('Name','Contraintes max sigxx 2D-Q4 (rouge) / poutre (bleu)')   
hold on,
% La contrainte sigxx poutre est calculee au centre des elements de la rangee du bas
x=0:L; y=Mf*(1-1/ney)*x/(2*I); plot(x,y,'b');
dx=L/(10*nex);
for i=1:nex      %contrainte sigxx sur les elements de la rangee du bas
    x1=(i-1)*L/nex ; x2 = i*L/nex ;  
    sign=sige(1+ney*(i-1));
    x=x1:dx:x2; y=sign*x/x; plot(x,y,'r'),
end
grid

