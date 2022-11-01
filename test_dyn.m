%% clear
clear all; close all; clc;
%% load model
path_to_urdf = 'urdf/humanoid.urdf';
joint_num = 5;

robot = importrobot(path_to_urdf);
robot.DataFormat = 'column';
robot.Gravity = [0, 0, -9.81];
%%
iter = 100;
for i = 1:iter
    q = -2*pi + 4*pi*rand(joint_num,1);
    %q = zeros(joint_num,1)+pi/10;
    q_d = zeros(joint_num,1)+pi/10;
    q_2d = zeros(joint_num,1)+pi/10;

    tau_matlab = inverseDynamics(robot,q,q_d,q_2d);
    tau_lgr = M_mat(q)*q_2d + ...
                C_mat(q, q_d)*q_d + ...
                G_mat(q);
            
    assert(norm(tau_matlab - tau_lgr) < 1e-8);
end
fprintf("Rigid Body Inverse Dynamics Test - OK!\n");
%%
% q = zeros(joint_num,1)+pi/10;
%% test M_mat, G_mat
% M_mat_matlab = massMatrix(robot,q);
% M_mat_lgr = M_mat(q);
% gravTorq = gravityTorque(robot,q);
% gravTorq_lgr = G_mat(q);
%% test trans mat
% source_id = 5;
% target_body_name = robot.Base.Name;
% source_body_name = robot.Bodies{source_id}.Name;
% transform_matlab = getTransform(robot, q, source_body_name, target_body_name);
% trans_all = Trans_mat(q);
% transform_matlab_lgr = trans_all(:,:,source_id);
%% test jacobian
% end_id = 5;
% end_name = robot.Bodies{end_id}.Name;
% jacobian_matlab = geometricJacobian(robot,q,end_name);
% jacobian_lgr_all = geometric_Jacobian(q);
% jacobian_lgr = jacobian_lgr_all(:,:,end_id);