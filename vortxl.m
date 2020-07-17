
% ·····················································
% Velocity induced by a constant strenght linear vortex
% ·····················································

function [u] = vortxl (X1,X2,XP,gamma)

x1 = X1(1) ; y1 = X1(2) ; z1 = X1(3) ;
x2 = X2(1) ; y2 = X2(2) ; z2 = X2(3) ;
xp = XP(1) ; yp = XP(2) ; zp = XP(3) ;


r0 = X2-X1 ;
r1 = XP-X1 ;
r2 = XP-X2 ;

norm_r1 = norm(r1) ;
norm_r2 = norm(r2) ;

r1xr2 = cross(r1,r2) ;
norm_r1xr2 = norm(r1xr2) ;

tol = 1.0e-6 ; u = zeros(3,1) ;

if ( norm_r1>tol && norm_r2>tol && norm_r1xr2>tol )
    
    inv_r1xr2 = 1.0 / norm_r1xr2 ;
    inv_r1 = 1.0 / norm_r1 ;
    inv_r2 = 1.0 / norm_r2 ;
    
    a = r0 * inv_r1xr2 ;
    b = r1*inv_r1 - r2*inv_r2 ;
    c = dot(a,b) ;
    
    u = gamma*0.25/pi*c*r1xr2*inv_r1xr2 ;

else
    u = 0.0 ;
end
    
    

% 
% 
% 
% 
% 
% 
% 
% tol = 1.0e-6 ; u = zeros(3,1) ;
% 
% a  =  (yp-y1)*(zp-z2) - (zp-z1)*(yp-y2) ;
% b  = -(xp-x1)*(zp-z2) + (zp-z1)*(xp-x2) ;
% c  =  (xp-x1)*(yp-y2) - (yp-y1)*(xp-x2) ;
% 
% d  = a*a + b*b + c*c ;
% r1 = sqrt( (xp-x1)*(xp-x1) + (yp-y1)*(yp-y1) + (zp-z1)*(zp-z1) ) ;
% r2 = sqrt( (xp-x2)*(xp-x2) + (yp-y2)*(yp-y2) + (zp-z2)*(zp-z2) ) ;
% 
% if ( d>tol & r2>tol & r1>tol ) 
%     ror1 = (x2-x1)*(xp-x1)+(y2-y1)*(yp-y1)+(z2-z1)*(zp-z1) ;
%     ror2 = (x2-x1)*(xp-x2)+(y2-y1)*(yp-y2)+(z2-z1)*(zp-z2) ;
%     com  = (gamma/(4.*pi*d))*((ror1/r1)-(ror2/r2)) ;
%     u(1) = a * com ;
%     u(2) = b * com ;
%     u(3) = c * com ;
% else
%     u(1) = 0.0 ;
%     u(2) = 0.0 ;
%     u(3) = 0.0 ;
% end

end