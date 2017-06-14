% Trilateration algorithm
% paper "An algebraic solution to the multilateration problem"
% Author: Norrdine, Abdelmoumen  (norrdine@hotmail.de)
% https://www.researchgate.net/publication/275027725_An_Algebraic_Solution_to_the_Multilateration_Problem
% usage: [N1 N2] = Trilateration(P,S,W) 
% P = [P1 P2 P3 P4 ..] Reference points matrix
% S = [s1 s2 s3 s4 ..] distance matrix.
% W : Weights Matrix (Statistics).
% N : calculated solution
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY!! 

function [N1, N2] = TrilaterationRub(refPoints, distances)
% refPoints. numPoints x 3 matrix.
% distances. numPoints-element column vector

numPoints = numel(distances);

A = [ones(numPoints, 1), -2*refPoints];
b = distances.^2 - sum(A.^2, 2);

if (numPoints == 3)
%     warning off;
    Xp=pinv(A)*b; 
    xp = Xp(2:4,:);
    Z = null(A,'r');
    z = Z(2:4,:);
    if rank (A)==3
        %Polynom coeff.
        a2 = z(1)^2 + z(2)^2 + z(3)^2 ;
        a1 = 2*(z(1)*xp(1) + z(2)*xp(2) + z(3)*xp(3))-Z(1);
        a0 = xp(1)^2 +  xp(2)^2+  xp(3)^2-Xp(1);
        p = [a2 a1 a0];
        t = roots(p);

        %Solutions
        N1 = (Xp + t(1)*Z)';
        N2 = (Xp + t(2)*Z)';
    end
end

if  (numPoints > 3)
    Xpdw=pinv(A)*b;
 
    N1 = Xpdw;
    N2 = N1;
end