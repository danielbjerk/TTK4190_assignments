syms p q r m L B H

r_bg = [-3.7; 0; H / 2];
nu_2 = [p; q; r];


I_bg = sym('I_bg', [3, 3]);
I_bg(:, :) = 0;
I_bg(3, 3) = (1/12) * m * (L^2 + B^2);

I_bb = I_bg - m * (skew(r_bg))^2;

MRB = [m*eye(3) m*skew(r_bg);
       m*skew(r_bg) I_bb];
   
MRB_DOF_126 = MRB;
MRB_DOF_126(3:5, :) = [];
MRB_DOF_126(:, 3:5) = [];
disp("MRB for DOF 1, 2 and 6: ");
disp(MRB_DOF_126);

CRB = [m * skew(nu_2), -m * skew(nu_2) * skew(r_bg);
       m * skew(r_bg) * skew(nu_2), -skew(I_bb * nu_2)];

CRB_DOF_126 = CRB;
CRB_DOF_126(3:5, :) = [];
CRB_DOF_126(:, 3:5) = [];
disp("CRB for DOF 1, 2 and 6: ");
disp(CRB_DOF_126);


function S = skew(x)
    S=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
end