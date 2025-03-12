
format long
format short

% To run this code, you need to have the SeDuMi solver installed.

% Solve Lasserre hierarchy

% Line 4 of Table I in the Ms:
%alpha12 = 54*2*pi/360;
%alpha13 = 112*2*pi/360;
%alpha23 = 194*2*pi/360;

% Line 3 of Table I in the Ms:
alpha12 = 58.4*2*pi/360;
alpha13 = 121.6*2*pi/360;
alpha23 = 180*2*pi/360;

 % Choose Task
   No_Of_Task=2;
    
    switch No_Of_Task
        case 1,
            % Task 1:
            % (|n1-n2| = |n1-n3|)
             % Q = Q^123_mirror
            om12 = (1/2)*(1+(1-abs(sin(alpha12)))/cos(alpha12));
            om13 = (1/2)*(1+(1-abs(sin(alpha13)))/cos(alpha13));
            om23 = (1/2)*(1+(1-abs(sin(alpha23)))/cos(alpha23));
        case 2,
            % Task 2: 
            % (|n1-n2| = |n2-n3|)
            % Q = Q^213_mirror
            om12 = (1/2)*(1+(1-abs(sin(alpha12)))/cos(alpha12));
            om23 = (1/2)*(1+(1-abs(sin(alpha13)))/cos(alpha13));
            om13 = (1/2)*(1+(1-abs(sin(alpha23)))/cos(alpha23));
        case 3,
            % Task 3: 
            % (|n1-n3| = |n2-n3|)
            % Q = Q^312_mirror
            om23 = (1/2)*(1+(1-abs(sin(alpha12)))/cos(alpha12));
            om13 = (1/2)*(1+(1-abs(sin(alpha13)))/cos(alpha13));
            om12 = (1/2)*(1+(1-abs(sin(alpha23)))/cos(alpha23));
    end

% define variables
n = sdpvar(3,3);
m = sdpvar(6,3);

% maximize Q
% subject to the following constraints
F = [];
for i=1:3
F = F + [n(i,:)*n(i,:)'<=1];
end;
for i=1:6
F = F + [m(i,:)*m(i,:)'==1];
end;
% The symmetry constraint defined by Eq. (30) in the Ms:
F = F + [n(2,:)*n(2,:)' - 2*n(1,:)*n(2,:)' ==  n(3,:)*n(3,:)' - 2*n(1,:)*n(3,:)'];

om = om12 + om13 + om23; 

% define Q (according to Eq.(29) in the Ms):
Q = om + om12*(n(1,:)+n(2,:))*m(1,:)' + (1-om12)*(n(1,:)-n(2,:))*m(2,:)' + ...
    om13*(n(1,:)+n(3,:))*m(3,:)' + (1-om13)*(n(1,:)-n(3,:))*m(4,:)' + ...
    om23*(n(2,:)+n(3,:))*m(5,:)' + (1-om23)*(n(2,:)-n(3,:))*m(6,:)';

obj = real(-Q);

ops= sdpsettings('solver','sedumi');

level = 1;
[sol,hij,momentdata]=solvemoment(F,obj,ops,level);

format short

sol

relaxobj = -relaxvalue(obj)

return;
