%% clear
clear all; close all; clc;
%% add path
% add generate path
cur = pwd;
generate_path = fullfile(cur,'auto_gen');
if ~exist(generate_path, 'dir')
    mkdir(generate_path);
end
addpath(genpath('./'));

%% load model
path_to_urdf = 'urdf/humanoid.urdf';
robot = xml2struct(path_to_urdf);
% set joint num
joint_num = 5;
g = 9.81;

%% create generilized coordiates
q_sym = sym('q%d',[joint_num,1],'real');
qd_sym = sym('qd%d',[joint_num,1],'real');
q2d_sym = sym('q2d%d',[joint_num,1],'real');

%% compute transe mat T 
% Ti represent from joint_i frame to origin frame, P_0 = T*P_i
T = cell(joint_num, 1);
% convert all the T to mat format
T_all_mat = sym(zeros(4, 4, joint_num));
for i = 1:joint_num
    % the offset T
    rpy_i_offset = sym(str2num(robot.robot.joint{i}.origin.Attributes.rpy));
    Rot_i_offset = euler2rot(rpy_i_offset);
    Pos_i_offset = str2num(robot.robot.joint{i}.origin.Attributes.xyz)';  
    Trans_i_offset = sym([Rot_i_offset, Pos_i_offset; zeros(1,3), 1]);
    % when joint i revolute qi angle
    joint_axis = str2num(robot.robot.joint{i}.axis.Attributes.xyz)';
    Rot_i = Rot(q_sym(i),sym(joint_axis));
    Pos_i = sym(zeros(3,1));
    Trans_i = [Rot_i, Pos_i; sym(zeros(1,3)),sym(1)];
    % trans for joint i to parent joint, T_i = T_off*T_rot
    T_i_final = Trans_i_offset*Trans_i;
    % T_0 = T_1 * ... * T_i-1 * T_i
    if i>1
        T{i} = T{i-1}*T_i_final;
    else
        T{i} = T_i_final;
    end
    T_all_mat(:,:,i) = T{i};
end
%% compute Rot_mat
R = cell(joint_num, 1);
for i = 1:joint_num
    R{i} = T{i}(1:3,1:3);
end
%% compute angular velocity jacobian Jw
Jw = cell(joint_num, 1);
for i = 1:joint_num
    temp = [];
    for j = 1:i
        joint_axis = str2num(robot.robot.joint{j}.axis.Attributes.xyz)';
        temp = [temp, R{j}*joint_axis];
    end
    Jw{i} = [temp, repmat([0;0;0],1,joint_num-i)];
end

%% compute linear velocity jacobian Jv for each link com
Jv = cell(joint_num, 1);
Com = cell(joint_num, 1);
for i = 1:joint_num
    com_pos_i = str2num(robot.robot.link{i+1}.inertial.origin.Attributes.xyz);
    Pos_i = T{i}*[com_pos_i,1]';
    %Pos_i = T{i}*[0,0,0,1]';
    Pos_i = Pos_i(1:3);
    Jv{i} = jacobian(Pos_i, q_sym);
    Com{i} = Pos_i;
end
%%
geometric_Jacobian = sym(zeros(6,joint_num,joint_num));
for i = 1:joint_num
    geometric_Jacobian(1:3,:,i) = Jw{i};
    geometric_Jacobian(4:6,:,i) = Jv{i};
end
%% get I_mat
I_mat = cell(joint_num, 1);
Link_mass = cell(joint_num, 1);
for i = 1:joint_num
    ixx = str2double(robot.robot.link{i+1}.inertial.inertia.Attributes.ixx);
    ixy = str2double(robot.robot.link{i+1}.inertial.inertia.Attributes.ixy);
    ixz = str2double(robot.robot.link{i+1}.inertial.inertia.Attributes.ixz);
    iyy = str2double(robot.robot.link{i+1}.inertial.inertia.Attributes.iyy);
    iyz = str2double(robot.robot.link{i+1}.inertial.inertia.Attributes.iyz);
    izz = str2double(robot.robot.link{i+1}.inertial.inertia.Attributes.izz);
    I_mat{i} = [ixx, ixy, ixz;
                ixy, iyy, iyz;
                ixz, iyz, izz]; 
    Link_mass{i} = str2double(robot.robot.link{i+1}.inertial.mass.Attributes.value);
end

%% get M_mat 
M_mat = zeros(joint_num, joint_num);
for i = 1:joint_num
    M_mat = M_mat + Link_mass{i}*Jv{i}'*Jv{i}+Jw{i}'*R{i}*I_mat{i}*R{i}'*Jw{i};
end

%% get G_mat
P = 0;
for i = 1:joint_num
    P = P+Link_mass{i}*g*Com{i}(3);
end
G_mat = jacobian(P, q_sym)';

%% get C_mat
c = sym(zeros(joint_num,joint_num,joint_num));
for k = 1:joint_num
    for i = 1:joint_num
        for j =1:joint_num
            c(i,j,k) = 0.5*(diff(M_mat(k,j),q_sym(i))+...
                            diff(M_mat(k,i),q_sym(j))-...
                            diff(M_mat(i,j),q_sym(k)));
        end
    end
end

C_mat = sym(zeros(joint_num, joint_num));
for k = 1:joint_num
    for j = 1:joint_num
        for i = 1:joint_num
            C_mat(k, j) = C_mat(k, j) + c(i,j,k)*qd_sym(i);
        end
    end
end
%% generate function
matlabFunction(M_mat, 'File','auto_gen/M_mat',...
               'Vars',{q_sym}, 'Optimize', true);  
%%
matlabFunction(C_mat, 'File','auto_gen/C_mat',...
               'Vars',{q_sym, qd_sym}, 'Optimize', false);
matlabFunction(G_mat, 'File','auto_gen/G_mat',...
               'Vars',{q_sym}, 'Optimize', true);
%% generate test function
matlabFunction(T_all_mat, 'File','auto_gen/Trans_mat',...
               'Vars',{q_sym}, 'Optimize', true);

%% 
matlabFunction(geometric_Jacobian, 'File','auto_gen/geometric_Jacobian',...
               'Vars',{q_sym}, 'Optimize', true);
