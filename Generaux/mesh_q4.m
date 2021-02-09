function [Coord,Connec]=mesh_q4(nex,ney,a,b)
%
% Maillage réglé d'un rectangle (a,b) en (nex,ney)éléments Q4
% H.Oudin
%  
Coord=[];
for j=0:ney
    for i=0:nex
        Coord=[Coord;[i*a/nex  j*b/ney]];
    end
end
            
r=0;
Nprop=[];
for n=1:nex*ney
    Connec(n , : ) = [ n+r  n+r+1  n+r+nex+2  n+r+nex+1];
    r=floor(n/nex);
end
end