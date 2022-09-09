function rot_3D = rotation_matrix_3D(chi, psi, theta)
if nargin == 1
     psi = chi(2); theta = chi(3);chi = chi(1);
end
Rx = [ 1,       0,      0 ;
       0, cos(chi), -sin(chi);
       0, sin(chi), cos(chi)];
  
Ry = [ cos(psi),   0,  sin(psi) ;
       0,     1,       0;
       -sin(psi), 0, cos(psi)];

Rz = [ cos(theta),   -sin(theta),  0 ;
       sin(theta), cos(theta), 0;
       0,            0,       1];
 
rot_3D = Rx*Ry*Rz;


end