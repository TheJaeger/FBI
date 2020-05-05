function aDT = axonal_DT(clm)

% re-define the fODF coefficients for ease of developing the
       % intra-axonal diffusion tensor, A (below)
                                            c00 = clm(1);
c2_2 = clm(2); c2_1 = clm(3); c20 = clm(4); c21 = clm(5); c22 = clm(6);

% Fist, the 6 unique elements of the diffusion tensor...
A11 = ((sqrt(30)/3)*c00 - (sqrt(6)/3)*c20 + c22 + c2_2); 
A22 = (sqrt(30)/3)*c00 - (sqrt(6)/3)*c20 - c22 - c2_2; 
A33 = ((sqrt(30)/3)*c00 + (2*sqrt(6)/3)*c20);
A12 = 1i*(c22 - c2_2); % = A21
A13 = (-c21 + c2_1); % = A31
A23 = 1i*(-c21 - c2_1); % = A32

% Note, DT_a is symmetric so only need 6 unique elements
scale_factor = 1./(sqrt(30)*repmat(c00,[1,6]));                 
aDT = scale_factor.*[A11 A22 A33 A12 A13 A23];                     
aDT = real(aDT);
aDT(isnan(aDT)) = 0; 
aDT = aDT'; 

end