clear all; %#ok<CLALL>
clc;
%% Variable Definitions
%World properties
syms g t a_vel E_kin E_trans E_rot E_pot fdyn;
%g=9.81;
%Link properties
link_types = ['r', 'r', 'r'];
links = size(link_types, 2);

%Link properties
m =[0.17729193;0.22664436;0.29046611]%sym('m', [1 links]);
c =[0.03836329;0.16404998;0.13002176]% sym('c', [1 links]);

%Angle functions
theta   = sym([0, 0, 0]);
dtheta  = sym([0, 0, 0]);
ddtheta = sym([0, 0, 0]);

for i=1:links
    if link_types(i) == 'r'
        theta(i) = str2sym(strcat('theta_', num2str(i), '(t)'));
        dtheta(i) = str2sym(strcat('dtheta_', num2str(i), '(t)'));
        ddtheta(i) = str2sym(strcat('ddtheta_', num2str(i), '(t)'));
    else
        theta(i) = str2sym(strcat('theta_', num2str(i)));
    end
end

L = [0.05616;0.21869;0.26964]%sym(zeros(links, 1));
dL = sym(zeros(links, 1));
ddL = sym(zeros(links, 1));
%for i=1:links
%    if link_types(i) == 'p'
%        L(i) = str2sym(strcat('L', num2str(i), '(t)'));
%        dL(i) = str2sym(strcat('L', num2str(i), '(t)'));
%        ddL(i) = str2sym(strcat('L', num2str(i), '(t)'));
%    else
%        L(i) = str2sym(strcat('L', num2str(i)));
%    end
%end

%Inertial matrices
%sym('I1', [3 3]);
I(:,:,1)=[
     0.00035613,       0,                           -0.00000003;
     0,                0.00033843,                   0.00000530;
     -0.00000003,      0.00000530,                   0.00004754;
 ];

I(:,:,2)=[
     0.00006272,               -0.00002635,                -0.00001513;
    -0.00002635,                0.00704305,                 0.00000022;
    -0.00001513,                0.00000022,                 0.00702286;
 ];

I(:,:,3)=[
     0.00016544,   0.00000588,                 0;
    0.00000588,    0.00572234,                -0.00000025;
     0,            -0.00000025,               0.00562428;
 ];

%Rotate inertial matrices to line up with kinematic coordinates
% IR = [
%     0,  1,  0;
%     1,  0,  0;
%     0,  0,  1;
% ];
% for i=1:links
%     I(:,:,i) = IR*I(:,:,i)*IR.';
% end
% 
% clearvars IR

%DH parameters 
DH = [
        0,       0,         theta(1),        L(1);
       -pi/2,       0,       theta(2)-pi/2,     0;
        0,       L(2),      theta(3),           0;
        0,       L(3),             0,           0;
];



%% Forward Kinematics
fkin = forward_kinematics(DH);

%% Forward Dynamics
%Define position functions for links and COM
x =  sym(zeros(3, 1));
xc = sym(zeros(3, 1));
v =  sym(zeros(3, 1));
vc = sym(zeros(3, 1));
for i=1:links
     x(:,:,i) = trans_matrix(fkin, i+1);
    xc(:,:,i) = subs(x(:,:,i), strcat('L',num2str(i)), c(i));
    
     v(:,:,i) = diff( x(:,:, i), t);
    vc(:,:,i) = diff(xc(:,:, i), t);
   
    for j=1:links
         v(:,:,i) = subs(  v(:,:,i), diff(theta(j), t), dtheta(j) );
        vc(:,:,i) = subs( vc(:,:,i), diff(theta(j), t), dtheta(j) );
    end
end

a_vel = 0;    E_trans = 0;    E_rot = 0;    E_pot = 0;
for i=1:links
    %Calculate new total angular velocity around base frame
    a_vel = a_vel + rot_matrix(fkin, i)*[0;0;1]*diff(theta(i), t);
    a_vel = simplify(a_vel, 20);

    %Calculate translational energy
    E_trans = E_trans + 0.5*(m(i)*(vc(:, :, i).'*vc(:, :, i)));
        %disp(0.5*(m(i)*(vc(:, :, i).'*vc(:, :, i))))
    E_trans = simplify(E_trans, 20);
    
    %Calculate rotational energy
    E_rot = E_rot + 0.5*(a_vel).'*I(:,:,i)*(a_vel);
        %disp(0.5*(a_vel).'*I(:,:,i)*(a_vel))
    E_rot = simplify(E_rot, 20);
    
    %Calculate potential energy
    E_pot = E_pot + g*m(i)*xc(3,1,i);
    E_pot = simplify(E_pot, 20);
end

lagrange = E_trans + E_rot - E_pot;

fdyn = sym('o', [links 1]);
for i=1:links
    fdyn(i) = diff(functionalDerivative(lagrange, dtheta(i)), t) - functionalDerivative(lagrange, theta(i));
    for j=1:links
         fdyn(i) = subs( fdyn(i),    diff(theta(j), t),  dtheta(j) );
         fdyn(i) = subs( fdyn(i), diff(theta(j), t, t), ddtheta(j) );
         fdyn(i) = subs( fdyn(i),   diff(dtheta(j), t), ddtheta(j) );
    end
    fdyn(i) = vpa(simplify(fdyn(i), 50),3);
end

for i=1:links
    disp(vpa(fdyn(i),5));
end