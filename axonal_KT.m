function aKT = axonal_KT(C_complex)
%
% Hunter Moss
% beta version
% 08/29/2017
%
% Following is the axonal Kurtosis (beta version 08/22/17 from Jens) calculated from FBI:
% 
% The aKT is defined as,
% 
%               W^(a) = W^(a)_ijkl = 27*U_ijkl - 9*(V_ij*V_kl - V_ik*V_jl - V_il*V_jk)  % Eq.[3] 
% Re-define the fODF coefficients for ease of developing the intra-axonal kurtosis tensor (aKT), W (below)
%
                                                                                            c00 = C_complex(:,1);
                                              c2_2 = C_complex(:,2); c2_1 = C_complex(:,3); c20 = C_complex(:,4); c21 = C_complex(:,5); c22 = C_complex(:,6);
c4_4 = C_complex(:,7); c4_3 = C_complex(:,8); c4_2 = C_complex(:,9); c4_1 = C_complex(:,10);c40 = C_complex(:,11);c41 = C_complex(:,12);c42 = C_complex(:,13); c43 = C_complex(:,14); c44 = C_complex(:,15);
%
% 
% We fist define explicitly all unique the elements of V,
% 
V11 = (1./(10*c00)).*((10/3)*c00 - (2*sqrt(5)/3)*c20 + (sqrt(30)/3)*(c22 + c2_2)); % Eq.[4]
V22 = (1./(10*c00)).*((10/3)*c00 - (2*sqrt(5)/3)*c20 - (sqrt(30)/3)*(c22 + c2_2)); % Eq.[5]
V33 = (1./(10*c00)).*((10/3)*c00 + (4*sqrt(5)/3)*c20);                            % Eq.[6]
V12 = (1i./(10*c00)).*((sqrt(30)/3)*(c22 - c2_2));                                % Eq.[7]
V13 = -(1./(10*c00)).*((sqrt(30)/3)*(c21 - c2_1));                                % Eq.[8]
V23 = -(1i./(10*c00)).*((sqrt(30)/3)*(c21 + c2_1));                               % Eq.[9]
% 
% and of U,
% 
U1111 = (1./(10*c00)).*(2*c00 - (4*sqrt(5)/7)*c20 + (2*sqrt(30)/7)*(c22 + c2_2) + (2/7)*c40 - (2*sqrt(10)/21)*(c42 + c4_2) + (sqrt(70)/21)*(c44 + c4_4)); % Eq.[10]
U2222 = (1./(10*c00)).*(2*c00 - (4*sqrt(5)/7)*c20 - (2*sqrt(30)/7)*(c22 + c2_2) + (2/7)*c40 + (2*sqrt(10)/21)*(c42 + c4_2) + (sqrt(70)/21)*(c44 + c4_4)); % Eq.[11]
U3333 = (1./(10*c00)).*(2*c00 + (8*sqrt(5)/7)*c20 + (16/21)*c40);                                                                                         % Eq.[12]
U1112 = (1i./(10*c00)).*((sqrt(30)/7)*(c22 - c2_2) - (sqrt(10)/21)*(c42 - c4_2) + (sqrt(70)/21)*(c44 - c4_4));                                            % Eq.[13]
U1113 = (1./(10*c00)).*(-(sqrt(30)/7)*(c21 - c2_1) + (sqrt(5)/7)*(c41 - c4_1) - (sqrt(35)/21)*(c43 - c4_3));                                              % Eq.[14]
U1222 = (1i./(10*c00)).*((sqrt(30)/7)*(c22 - c2_2) - (sqrt(10)/21)*(c42 - c4_2) - (sqrt(70)/21)*(c44 - c4_4));                                            % Eq.[15]
U2223 = -(1i./(10*c00)).*((sqrt(30)/7)*(c21 + c2_1) - (sqrt(5)/7)*(c41 + c4_1) - (sqrt(35)/21)*(c43 + c4_3));                                             % Eq.[16]
U1333 = (1./(10*c00)).*(-(sqrt(30)/7)*(c21 - c2_1) - (4*sqrt(5)/21)*(c41 - c4_1));                                                                        % Eq.[17]
U2333 = -(1i./(10*c00)).*((sqrt(30)/7)*(c21 + c2_1) + (4*sqrt(5)/21)*(c41 + c4_1));                                                                       % Eq.[18]
U1122 = (1./(10*c00)).*((2/3)*c00 - (4*sqrt(5)/21)*c20 + (2/21)*c40 - (sqrt(70)/21)*(c44 + c4_4));                                                        % Eq.[19]
U1133 = (1./(10*c00)).*((2/3)*c00 + (2*sqrt(5)/21)*c20 + (sqrt(30)/21)*(c22 + c2_2) - (8/21)*c40 + (2*sqrt(10)/21)*(c42 + c4_2));                         % Eq.[20]
U2233 = (1./(10*c00)).*((2/3)*c00 + (2*sqrt(5)/21)*c20 - (sqrt(30)/21)*(c22 + c2_2) - (8/21)*c40 - (2*sqrt(10)/21)*(c42 + c4_2));                         % Eq.[21]
U1123 = -(1i./(10*c00)).*((sqrt(30)/21)*(c21 + c2_1) - (sqrt(5)/21)*(c41 + c4_1) + (sqrt(35)/21)*(c43 + c4_3));                                           % Eq.[22]
U1223 = (1./(10*c00)).*(-(sqrt(30)/21)*(c21 - c2_1) + (sqrt(5)/21)*(c41 - c4_1) + (sqrt(35)/21)*(c43 - c4_3));                                            % Eq.[23]
U1233 = (1i./(10*c00)).*((sqrt(30)/21)*(c22 - c2_2) + (2*sqrt(10)/21)*(c42 - c4_2));                                                                      % Eq.[24]

% Now, the 15-independent elements of the intra-axonal KT may be built with combinations of the above formulas:
W1111 = 27*U1111 - 27*(V11.*V11);
W2222 = 27*U2222 - 27*(V22.*V22);
W3333 = 27*U3333 - 27*(V33.*V33);

W1112 = 27*U1112 - 18*(V11.*V12) - 9*(V12.*V11);
W1113 = 27*U1113 - 18*(V11.*V13) - 9*(V13.*V11);
W1222 = 27*U1222 - 27*(V12.*V22);

W2223 = 27*U2223 - 18*(V22.*V23) - 9*(V23.*V22);
W1333 = 27*U1333 - 27*(V13.*V33);
W2333 = 27*U2333 - 27*(V23.*V33);

W1122 = 27*U1122 - 18*(V12.*V12) - 9*(V11.*V22);
W1133 = 27*U1133 - 18*(V13.*V13) - 9*(V11.*V33);
W2233 = 27*U2233 - 18*(V23.*V23) - 9*(V22.*V33);
W1123 = 27*U1123 - 9*(V11.*V23)- 9*(V12.*V13) - 9*(V13.*V12);
W1223 = 27*U1223 - 18*(V12.*V23) - 9*(V13.*V22); 
W1233 = 27*U1233 - 18*(V12.*V23) - 9*(V12.*V33);

% Combine these 15-elements into a vector
aKT = [W1111,W2222,W3333,W1112,W1113,W1222,W1333,W2223,W2333,W1122,W1133,W2233,W1123,W1223,W1233];
aKT = real(aKT);
aKT(isnan(aKT)) = 0; 
aKT = aKT';

end


