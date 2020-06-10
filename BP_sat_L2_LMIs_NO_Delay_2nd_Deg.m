clc
clear

% %
% % figure;
% % set(gca,'Linewidth',2);
% % plot(u(:,1),u(:,2),'Linewidth',2);
% % set(0,'DefaultTextInterpreter','Latex');set(gca,'FontSize',26);232
% % set(gca,'TicklabelInterpreter','Latex')
% % set( get( gca, 'XLabel' ), 'Interpreter', 'Latex' );
% % set( get( gca, 'YLabel' ), 'Interpreter', 'Latex' );
% % title('Control effort','FontSize', 26);
% % xlabel('Time [sec]','FontSize', 26);
% % ylabel('u(t)','FontSize', 26);

addpath(genpath('/home/saeed/Downloads/YALMIP'));
addpath(genpath('/home/saeed/Downloads/mosek/9.0/toolbox/r2015a/'));
%run /home/saeed/Downloads/cvx/cvx_startup.m
yalmip('clear')
% load('EngSpeed.mat')
% load('GridData.mat')

inequality=1e-8;
kappa=5;%4.9;%
delta_1=1e-12; % Diagonal small term

W_e=.9;%.2; Tracking error weight
W_u=.01;%2;% Input weight

[A_p,B_p,D_p,C_yp,C_zp,D_zu,D_zw]=systemdynamics(1,1,1,1,1);%systemdynamics(rho_1,rho_2,delta_1,W_e,W_u);
n=size(A_p,1);
n_w=size(D_p,2);
n_u=size(B_p,2);
n_z=size(C_zp,1);
n_y=size(C_yp,1);

%%% Saturation
u_bar=120;
delta=4e-5;% snorm(w(t),2)^-2
control=1;

%%%% Uncertainty
epsilon=300;%289;%
E_A=[.02,0;
    0,0];
G_A=[.005,0;
    0,0];

tau_min=100; %%% rho_1 Lag.
tau_max=200;
v_1=.06;

K_min = .45; %%% rho_2 Sens.
K_max=.66;
v_2=2e-4;

theta_min = 0.01; %%% rho_3 delay
theta_max=8;
v_3=.5;
theta_prime=1;

% step_rho_1=(tau_max-tau_min)/4;
% step_rho_2=(K_max-K_min)/4;
% step_rho_3=(theta_max-theta_min)/4;

for jj=1:size(epsilon,2)
    for j=1:size(kappa,2)
       
        P_bar_0=sdpvar(2*n);
        P_bar_1=sdpvar(2*n);
        P_bar_2=sdpvar(2*n);
        P_bar_3=sdpvar(2*n);
        P_bar_4=sdpvar(2*n);
      
        Q_bar_0=sdpvar(2*n);
        Q_bar_1=sdpvar(2*n);
        Q_bar_2=sdpvar(2*n);
        Q_bar_3=sdpvar(2*n);
        Q_bar_4=sdpvar(2*n);
       
        S_bar_0=sdpvar(2*n);
        S_bar_1=sdpvar(2*n);
        S_bar_2=sdpvar(2*n);
        S_bar_3=sdpvar(2*n);
        S_bar_4=sdpvar(2*n);
        
        W_bar=sdpvar(2*n);
        S1_bar=sdpvar(2*n,2*n,'full');
        gamma_sq=sdpvar(1);
        
        beta=sdpvar(1);
        
        T_bar_0= diag(sdpvar(n_u,1));
        T_bar_1= diag(sdpvar(n_u,1));
        T_bar_2= diag(sdpvar(n_u,1));
        T_bar_3= diag(sdpvar(n_u,1));
        T_bar_4= diag(sdpvar(n_u,1));
        
        T_1_bar_0= diag(sdpvar(n_u,1));
        T_1_bar_1= diag(sdpvar(n_u,1));
        T_1_bar_2= diag(sdpvar(n_u,1));
        T_1_bar_3= diag(sdpvar(n_u,1));
        T_1_bar_4= diag(sdpvar(n_u,1));
        
        X_0=sdpvar(n);
        X_1=sdpvar(n);
        X_2=sdpvar(n);
        X_3=sdpvar(n);
        X_4=sdpvar(n);
        
        Y_0=sdpvar(n);
        Y_1=sdpvar(n);
        Y_2=sdpvar(n);
        Y_3=sdpvar(n);
        Y_4=sdpvar(n);
        
        A_tilde_0=sdpvar(n,n,'full');
        A_tilde_1=sdpvar(n,n,'full');
        A_tilde_2=sdpvar(n,n,'full');
        A_tilde_3=sdpvar(n,n,'full');
        A_tilde_4=sdpvar(n,n,'full');
        
        B_tilde_0=sdpvar(n,n_y,'full');
        B_tilde_1=sdpvar(n,n_y,'full');
        B_tilde_2=sdpvar(n,n_y,'full');
        B_tilde_3=sdpvar(n,n_y,'full');
        B_tilde_4=sdpvar(n,n_y,'full');
       
        C_tilde_0=sdpvar(n_u,n,'full');
        C_tilde_1=sdpvar(n_u,n,'full');
        C_tilde_2=sdpvar(n_u,n,'full');
        C_tilde_3=sdpvar(n_u,n,'full');
        C_tilde_4=sdpvar(n_u,n,'full');
        
        D_K_0=sdpvar(n_u,n_y,'full');
        D_K_1=sdpvar(n_u,n_y,'full');
        D_K_2=sdpvar(n_u,n_y,'full');
        D_K_3=sdpvar(n_u,n_y,'full');
        D_K_4=sdpvar(n_u,n_y,'full');
        
        E_tilde_0=sdpvar(n,n_u,'full');
        E_tilde_1=sdpvar(n,n_u,'full');
        E_tilde_2=sdpvar(n,n_u,'full');
        E_tilde_3=sdpvar(n,n_u,'full');
        E_tilde_4=sdpvar(n,n_u,'full');
        
        A_tilde_theta_0=sdpvar(n,n,'full');
        A_tilde_theta_1=sdpvar(n,n,'full');
        A_tilde_theta_2=sdpvar(n,n,'full');
        A_tilde_theta_3=sdpvar(n,n,'full');
        A_tilde_theta_4=sdpvar(n,n,'full');
        A_tilde_theta_5=sdpvar(n,n,'full');
        A_tilde_theta_6=sdpvar(n,n,'full');
        
        B_tilde_theta_0=sdpvar(n,n_y,'full');
        B_tilde_theta_1=sdpvar(n,n_y,'full');
        B_tilde_theta_2=sdpvar(n,n_y,'full');
        B_tilde_theta_3=sdpvar(n,n_y,'full');
        B_tilde_theta_4=sdpvar(n,n_y,'full');
        B_tilde_theta_5=sdpvar(n,n_y,'full');
        B_tilde_theta_6=sdpvar(n,n_y,'full');
        
        E_tilde_theta_0=sdpvar(n,n_u,'full');
        E_tilde_theta_1=sdpvar(n,n_u,'full');
        E_tilde_theta_2=sdpvar(n,n_u,'full');
        E_tilde_theta_3=sdpvar(n,n_u,'full');
        E_tilde_theta_4=sdpvar(n,n_u,'full');
        E_tilde_theta_5=sdpvar(n,n_u,'full');
        E_tilde_theta_6=sdpvar(n,n_u,'full');
        
        G_tilde_0=sdpvar(n_u,2*n,'full');
        G_tilde_1=sdpvar(n_u,2*n,'full');
        G_tilde_2=sdpvar(n_u,2*n,'full');
        G_tilde_3=sdpvar(n_u,2*n,'full');
        G_tilde_4=sdpvar(n_u,2*n,'full');

        G_1_tilde_0=sdpvar(n_u,2*n,'full');
        G_1_tilde_1=sdpvar(n_u,2*n,'full');
        G_1_tilde_2=sdpvar(n_u,2*n,'full');
        G_1_tilde_3=sdpvar(n_u,2*n,'full');
        G_1_tilde_4=sdpvar(n_u,2*n,'full');
        
        LMI= blkdiag(W_bar,gamma_sq)>=inequality*eye(size(blkdiag(W_bar,gamma_sq)));
        LMI= LMI +(blkdiag([W_bar, S1_bar; S1_bar', W_bar],delta-beta)>=zeros(size(blkdiag([W_bar, S1_bar; S1_bar', W_bar],delta-beta))));%
        
        for rho_1 = 100:50:200%tau_min:step_rho_1:tau_max
            for rho_2 = [.45 .55 .66]%K_min:step_rho_2:K_max
                %                 for rho_3 = rho_3_min:step_rho_3:rho_3_max
                [A_p,B_p,D_p,C_yp,C_zp,D_zu,D_zw]=systemdynamics(rho_1,rho_2,delta_1,W_e,W_u);
                C_yp_delayed=C_yp;
                
                X= X_0 + rho_1 * X_1 + 1/2 * rho_1^2 * X_2 + rho_2 * X_3 + 1/2 * rho_2^2 * X_4;
                Y= Y_0 + rho_1 * Y_1 + 1/2 * rho_1^2 * Y_2 + rho_2 * Y_3 + 1/2 * rho_2^2 * Y_4;
                %                 X_1=0;
                %                 X_2=0;
                %                 Y_1=0;
                %                 Y_2=0;
                
                A_tilde = A_tilde_0 + rho_1 * A_tilde_1 + 1/2* rho_1^2 * A_tilde_2 + rho_2 * A_tilde_3 + 1/2* rho_2^2 * A_tilde_4;
                B_tilde = B_tilde_0 + rho_1 * B_tilde_1 + 1/2* rho_1^2 * B_tilde_2 + rho_2 * B_tilde_3 + 1/2* rho_2^2 * B_tilde_4;
                C_tilde = C_tilde_0 + rho_1 * C_tilde_1 + 1/2* rho_1^2 * C_tilde_2 + rho_2 * C_tilde_3 + 1/2* rho_2^2 * C_tilde_4;
                D_K = D_K_0 + rho_1 * D_K_1 + 1/2* rho_1^2 * D_K_2 + rho_2 * D_K_3 + 1/2* rho_2^2 * D_K_4;
                
                E_tilde = E_tilde_0 + rho_1* E_tilde_1 + 1/2* rho_1^2 * E_tilde_2 + rho_2 * E_tilde_3 + 1/2 * rho_2^2 * E_tilde_4;
                
                P_bar = P_bar_0 + rho_1 * P_bar_1 + 1/2 * rho_1^2 * P_bar_2 + rho_2 * P_bar_3 + 1/2 * rho_2^2 * P_bar_4;
                
                Q_bar = Q_bar_0 + rho_1 * Q_bar_1 + 1/2 * rho_1^2 * Q_bar_2 + rho_2 * Q_bar_3 + 1/2 * rho_2^2 * Q_bar_4;
                
                S_bar = S_bar_0 + rho_1 * S_bar_1 + 1/2 * rho_1^2 * S_bar_2 + rho_2 * S_bar_3 + 1/2 * rho_2^2 * S_bar_4;
                
                T_bar = T_bar_0 + rho_1 * T_bar_1 + 1/2 * rho_1^2 * T_bar_2 + rho_2 * T_bar_3 + 1/2 * rho_2^2 * T_bar_4;
                
                G_tilde = G_tilde_0 + rho_1 * G_tilde_1 + 1/2 * rho_1^2 * G_tilde_2 + rho_2 * G_tilde_3 + 1/2 * rho_2^2 * G_tilde_4;
                
                % beta =        beta_0 + rho * beta_1 + 1/2 * rho^2 * beta_2;
                
                for rho_1_d=100:50:200%tau_min:step_rho_1:tau_max
                    for rho_2_d = [.45 .55 .66]%K_min:step_rho_2:K_max
                        %                             for rho_3_d = rho_3_min:step_rho_3:rho_3_max
                        
                        A_tilde_theta = A_tilde_theta_0 + rho_1 * A_tilde_theta_1 + rho_2 * A_tilde_theta_2 +  rho_1 * rho_2 * A_tilde_theta_3 ...
                            + rho_1_d * A_tilde_theta_4 + rho_2_d * A_tilde_theta_5 + rho_1_d * rho_2_d * A_tilde_theta_6;
                        
                        B_tilde_theta = B_tilde_theta_0 + rho_1 * B_tilde_theta_1 + rho_2 * B_tilde_theta_2 + rho_1 * rho_2 * B_tilde_theta_3...
                            + rho_1_d * B_tilde_theta_4 + rho_2_d * B_tilde_theta_5 + rho_1_d * rho_2_d * B_tilde_theta_6;
                        
                        E_tilde_theta = E_tilde_theta_0 + rho_1 * E_tilde_theta_1 + rho_2 * E_tilde_theta_2 + rho_1 * rho_2 * E_tilde_theta_3 ...
                            + rho_1_d * E_tilde_theta_4 + rho_2_d * E_tilde_theta_5 + rho_1_d * rho_2_d * E_tilde_theta_6;
                        
                        C_tilde_delayed = C_tilde_0 + rho_1_d * C_tilde_1 + 1/2 * rho_1_d^2 * C_tilde_2 + rho_2_d * C_tilde_3 + 1/2 * rho_2_d^2 * C_tilde_4;
                        D_K_delayed = D_K_0 + rho_1_d * D_K_1 + 1/2* rho_1_d^2 * D_K_2 + rho_2_d * D_K_3 + 1/2* rho_2_d^2 * D_K_4;
                        
                        Q_bar_delayed = Q_bar_0 + rho_1_d * Q_bar_1 + 1/2 * rho_1_d^2 * Q_bar_2 + rho_2_d * Q_bar_3 + 1/2 * rho_2_d^2 * Q_bar_4;
                        T_1_bar_delayed = T_1_bar_0 + rho_1_d * T_1_bar_1 + 1/2 * rho_1_d^2 * T_1_bar_2 + rho_2_d * T_1_bar_3 + 1/2 * rho_2_d^2 * T_1_bar_4 ;
                        
                        G_1_tilde_delayed = G_1_tilde_0 + rho_1_d * G_1_tilde_1 + 1/2 * rho_1_d^2 * G_1_tilde_2 + rho_2_d * G_1_tilde_3 + 1/2 * rho_2_d^2 * G_1_tilde_4;
                        
                        P_bar_delayed = P_bar_0 + rho_1_d * P_bar_1 + 1/2 * rho_1_d^2 * P_bar_2 + rho_2_d * P_bar_3 + 1/2 * rho_2_d^2 * P_bar_4;
                        %         for counter=1:n_u
                        %             LMI=LMI+([beta, [C_tilde_delayed(counter,:), D_K_delayed(counter,:) * C_yp]-G_hat_delayed(counter,:);...
                        %                                 ([C_tilde(counter,:),D_K(counter,:) * C_yp]-G_hat(counter,:))', u_bar(counter)^2*P_bar_delayed]...
                        %                 >=zeros(size([beta, [C_tilde_delayed(counter,:), D_K_delayed(counter,:) * C_yp]-G_hat_delayed(counter,:);...
                        %                                 ([C_tilde(counter,:),D_K(counter,:) * C_yp]-G_hat(counter,:))', u_bar(counter)^2*P_bar_delayed])));
                        %         end
                        
                        %             for rho_h=800:step_h:3500
                        %                 S_bar_h = S_bar_0 + rho_h * S_bar_1 + 1/2 * rho_h^2 * S_bar_2;%%
                        S_bar_h=S_bar;
                        
                        A_hat_cl = [A_p*Y,A_p;A_tilde,X*A_p+B_tilde*C_yp];
                        A_hat_h = [B_p*C_tilde_delayed,B_p*D_K_delayed*C_yp_delayed;A_tilde_theta,B_tilde_theta*C_yp_delayed];
                        C_hat_h = D_zu*[C_tilde_delayed, D_K_delayed*C_yp_delayed];
                                    
                        for v_1=[-v_1,v_1]
                            for v_2=[-v_2,v_2]
                                for v_3=[-v_3,v_3]
                                    
                                    P_dot_bar = v_1*(P_bar_1 + rho_1*P_bar_2)+v_2*(P_bar_3 + rho_2*P_bar_4);
                                    
                                    L_11 = A_hat_cl + A_hat_cl'+ Q_bar + S_bar + P_dot_bar - W_bar;
                                    L_12 = P_bar - [Y,eye(n);eye(n),X] + kappa(j)*A_hat_cl';
                                    L_13 = S1_bar;
                                    L_14 = W_bar - S1_bar + A_hat_h;
                                    L_15 = [zeros(n,n_u);E_tilde] + G_tilde';
                                    L_16 = [-B_p * T_1_bar_delayed;E_tilde_theta];
                                    L_17 = [D_p;X*D_p];
                                    L_18 = [C_zp*Y,C_zp]';
                                    L_19 = [E_A,zeros(size(E_A,1),n);X*E_A,zeros(size(X,1),n)];
                                    L_110 = epsilon(jj)*[Y*G_A', zeros(size(Y,1),n);
                                        G_A' , zeros(size(G_A',1),n)];
                                    
                                    L_22 = -2*kappa(j)*[Y,eye(n);eye(n),X] + theta_max^2*W_bar;
                                    L_23 = zeros(2*n,2*n);
                                    L_24 = kappa(j)*A_hat_h;
                                    L_25 = kappa(j)*[zeros(n,n_u);E_tilde];
                                    L_26 = kappa(j)*[-B_p * T_1_bar_delayed;E_tilde_theta];
                                    L_27 = kappa(j)*[D_p;X*D_p];
                                    L_28 = zeros(2*n,size(C_zp',2));
                                    L_29 = kappa(j)*[E_A,zeros(size(E_A,1),n);X*E_A,zeros(size(X,1),n)];
                                    L_210 = zeros(2*n,2*n);
                                    
                                    L_33 = -S_bar_h - W_bar;
                                    L_34 = W_bar - S1_bar';
                                    L_35 = zeros(2*n,size(E_tilde,2));
                                    L_36 = zeros(2*n,size(T_1_bar_delayed,2));
                                    L_37 = zeros(2*n,size(D_p,2));
                                    L_38 = zeros(2*n,size(C_zp',2));
                                    L_39 = zeros(2*n,2*n);
                                    L_310 = zeros(2*n,2*n);
                                    
                                    
                                    L_44 = -2*W_bar+S1_bar+ S1_bar'-(1-v_3*theta_prime)*Q_bar_delayed;
                                    L_45 = zeros(2*n,size(E_tilde,2));
                                    L_46 = G_1_tilde_delayed';
                                    L_47 = zeros(2*n,size(D_p,2));
                                    L_48 = C_hat_h';
                                    L_49 = zeros(2*n,2*n);
                                    L_410 = zeros(2*n,2*n);
                                    
                                    
                                    L_55 = -2*T_bar;
                                    L_56 = zeros(size(T_bar,1),size(T_1_bar_delayed,2));
                                    L_57 = zeros(size(T_bar,1),size(D_p,2));
                                    L_58 = zeros(size(T_bar,1),size(C_zp',2));
                                    L_59 = zeros(size(T_bar,1),2*n);
                                    L_510 = zeros(size(T_bar,1),2*n);
                                    
                                    L_66= -2*T_1_bar_delayed;
                                    L_67= zeros(size(T_1_bar_delayed,1),size(D_p,2));
                                    L_68= -T_1_bar_delayed*D_zu';
                                    L_69= zeros(size(T_1_bar_delayed,1),2*n);
                                    L_610= zeros(size(T_1_bar_delayed,1),2*n);
                                    
                                    L_77= -eye(size(D_p,2));
                                    L_78= D_zw';
                                    L_79= zeros(size(D_zw',1),2*n);
                                    L_710= zeros(size(D_zw',1),2*n);
                                    
                                    L_88= -gamma_sq*eye(size(D_zw',2));
                                    L_89= zeros(size(L_88,1),2*n);
                                    L_810= zeros(size(L_88,1),2*n);
                                    
                                    L_99= -epsilon(jj)*eye(2*n);
                                    L_910= zeros(size(L_99,1),2*n);
                                    
                                    L_1010= -epsilon(jj)*eye(2*n);
                                    
                                    %                 L =[L_11,  L_12,  L_13,  L_14,  L_15,  L_16,   L_17;
                                    %                     L_12', L_22,  L_23,  L_24,  L_25,  L_26,   L_27;
                                    %                     L_13', L_23', L_33,  L_34,  L_35,  L_36,   L_37;
                                    %                     L_14', L_24', L_34', L_44,  L_45,  L_46,   L_47;
                                    %                     L_15', L_25', L_35', L_45', L_55,  L_56,   L_57;
                                    %                     L_16', L_26', L_36', L_46', L_56', L_66,   L_67;
                                    %                     L_17', L_27', L_37', L_47', L_57', L_67',  L_77];
                                    
                                    L =[L_11,   L_12,   L_13,   L_14,   L_15,   L_16,   L_17,   L_18,   L_19,   L_110;
                                        L_12',  L_22,   L_23,   L_24,   L_25,   L_26,   L_27,   L_28,   L_29,   L_210;
                                        L_13',  L_23',  L_33,   L_34,   L_35,   L_36,   L_37,   L_38,   L_39,   L_310;
                                        L_14',  L_24',  L_34',  L_44,   L_45,   L_46,   L_47,   L_48,   L_49,   L_410;
                                        L_15',  L_25',  L_35',  L_45',  L_55,   L_56,   L_57,   L_58,   L_59,   L_510;
                                        L_16',  L_26',  L_36',  L_46',  L_56',  L_66,   L_67,   L_68,   L_69,   L_610;
                                        L_17',  L_27',  L_37',  L_47',  L_57',  L_67',  L_77,   L_78,   L_79,   L_710;
                                        L_18',  L_28',  L_38',  L_48',  L_58',  L_68',  L_78',  L_88,   L_89,   L_810;
                                        L_19',  L_29',  L_39',  L_49',  L_59',  L_69',  L_79',  L_89',  L_99,   L_910;
                                        L_110', L_210', L_310', L_410', L_510', L_610', L_710', L_810', L_910', L_1010];
                                    
                                    
                                    
                                    %                 LMI=LMI+(blkdiag(P_bar,Q_bar,Q_bar_delayed,S_bar,S_bar_h,T_bar_delayed)>=...
                                    %                     1e-8*eye(size(blkdiag(P_bar,Q_bar,Q_bar_delayed,S_bar,S_bar_h,T_bar_delayed))));
                                    LMI=LMI+(L<=-inequality * eye(size(L)));
                                end
                                %                  LMI=LMI+(S_bar_h>=1e-8*eye(size(S_bar_h)));
                            end
                            %                                     LMI=LMI+(blkdiag(Q_bar_delayed,T_bar_delayed)>=1e-8*eye(size(blkdiag(Q_bar_delayed,T_bar_delayed))));
                        end
                        for counter=1:n_u
                            LMI_sat_d=[beta, [C_tilde_delayed(counter,:), D_K_delayed(counter,:) * C_yp]-G_1_tilde_delayed(counter,:);...
                                ([C_tilde_delayed(counter,:),D_K_delayed(counter,:) * C_yp]-G_1_tilde_delayed(counter,:))', u_bar(counter)^2*P_bar_delayed];
                            LMI=LMI + (LMI_sat_d>=0);
                        end
                        %                             end
                    end
                end
                for counter=1:n_u
                    LMI_sat=[beta, [C_tilde(counter,:), D_K(counter,:) * C_yp]-G_tilde(counter,:);...
                        ([C_tilde(counter,:),D_K(counter,:) * C_yp]-G_tilde(counter,:))', u_bar(counter)^2*P_bar];
                    LMI=LMI+(LMI_sat>=0);
                end
                LMI=LMI+(blkdiag(Q_bar,S_bar)>=inequality*eye(size(blkdiag(Q_bar,S_bar))));
                %                 end
            end
        end
        
        
        ops=sdpsettings('solver','mosek','debug',1);
%         ops.mosek.MSK_DPAR_PRESOLVE_TOL_X=1e-5;
%         ops.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-5;
%         ops.mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-5;
%         ops.mosek.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL=1e4/2;
        ops.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-5;
%         ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-6;  
        sol=optimize(LMI,gamma_sq,ops);
%         kappa_counter=j;
        kapp=kappa(j)
%         epsilon_counter=jj;
        epsil=epsilon(jj)
        sol
        
        gamma=sqrt(value(gamma_sq))
        if sol.problem==0
            AA(j,1,jj)=gamma;
            AA(j,2,jj)=kappa(j);
            AA(j,3,jj)=epsilon(jj);
            AA(j,4,jj)=sol.problem;
        else
            AA(j,1,jj)=1e5;
            AA(j,2,jj)=kappa(j);
            AA(j,3,jj)=epsilon(jj);
            AA(j,4,jj)=sol.problem;
        end
        if control==1
                T_min = 10;
                T_max = 35;
                K0 = 0.55;  % K0 = K_min+(K_max-K_min)*rand;
                K1 = 0.0045; %K1 = 0.002+0.005*rand;
                aK1 = 550;   % aK1 = 500+100*rand;
                T0 = 20;  % T0 = T_min+(T_max-T_min)*rand;
                aT2 = 10;  % 5+10*rand;
                aT1 = 10;   % aT1 = 5+10*rand;
                bT1 = 100;  % bT1 = 80+40*rand;
                tau0 = 150;
                tau1 = 5e-4;  % tau1 = 1e-3*rand;            
            S1_bar=value(S1_bar);S_bar=value(S_bar);W_bar=value(W_bar); P_bar=value(P_bar);Q_bar=value(Q_bar);
            X_0=value(X_0);X_1=value(X_1);X_2=value(X_2);X_3=value(X_3);X_4=value(X_4);
            
            Y_0=value(Y_0);Y_1=value(Y_1);Y_2=value(Y_2);Y_3=value(Y_3);Y_4=value(Y_4);
            
            A_tilde_0=value(A_tilde_0);A_tilde_1=value(A_tilde_1);A_tilde_2=value(A_tilde_2);A_tilde_3=value(A_tilde_3); A_tilde_4=value(A_tilde_4);
            
            A_tilde_theta_0=value(A_tilde_theta_0);A_tilde_theta_1=value(A_tilde_theta_1);A_tilde_theta_2=value(A_tilde_theta_2);
            A_tilde_theta_3=value(A_tilde_theta_3);A_tilde_theta_4=value(A_tilde_theta_4);A_tilde_theta_5=value(A_tilde_theta_5);
            A_tilde_theta_6=value(A_tilde_theta_6);
                  
            B_tilde_0=value(B_tilde_0);B_tilde_1=value(B_tilde_1);B_tilde_2=value(B_tilde_2);B_tilde_3=value(B_tilde_3);B_tilde_4=value(B_tilde_4);
            
            B_tilde_theta_0=value(B_tilde_theta_0);B_tilde_theta_1=value(B_tilde_theta_1);B_tilde_theta_2=value(B_tilde_theta_2);
            B_tilde_theta_3=value(B_tilde_theta_3);B_tilde_theta_4=value(B_tilde_theta_4);B_tilde_theta_5=value(B_tilde_theta_5);
            B_tilde_theta_6=value(B_tilde_theta_6);
            
            C_tilde_0=value(C_tilde_0);C_tilde_1=value(C_tilde_1);C_tilde_2=value(C_tilde_2);C_tilde_3=value(C_tilde_3);C_tilde_4=value(C_tilde_4);
            
            D_K_0=value(D_K_0);D_K_1=value(D_K_1);D_K_2=value(D_K_2);D_K_3=value(D_K_3);D_K_4=value(D_K_4);
            
            E_tilde_0=value(E_tilde_0);E_tilde_1=value(E_tilde_1);E_tilde_2=value(E_tilde_2);E_tilde_3=value(E_tilde_3);E_tilde_4=value(E_tilde_4);
            
            E_tilde_theta_0=value(E_tilde_theta_0);E_tilde_theta_1=value(E_tilde_theta_1);E_tilde_theta_2=value(E_tilde_theta_2);
            E_tilde_theta_3=value(E_tilde_theta_3);E_tilde_theta_4=value(E_tilde_theta_4);E_tilde_theta_5=value(E_tilde_theta_5);
            E_tilde_theta_6=value(E_tilde_theta_6);
            
            T_bar_0= value(T_bar_0);T_bar_1= value(T_bar_1);T_bar_2= value(T_bar_2);T_bar_3= value(T_bar_3);T_bar_4= value(T_bar_4);
            
            T_1_bar_0= value(T_1_bar_0);T_1_bar_1= value(T_1_bar_1);T_1_bar_2= value(T_1_bar_2);T_1_bar_3= value(T_1_bar_3);T_1_bar_4= value(T_1_bar_4);
        end
        beta=value(beta) % beta_0=value(beta_0);beta_1=value(beta_1);beta_2=value(beta_2);
        yalmip('clear')
    end
    
end

function [A_p,B_p,D_p,C_yp,C_zp,D_zu,D_zw]=systemdynamics(rho_1,rho_2,delta_1,W_e,W_u)
A_p=[-1/rho_1,0;
    -1,-delta_1];
B_p=[rho_2/rho_1;0];
D_p=[0,rho_2/rho_1;
    1,0];% Omega=[r;d]
%   D_p=[0,1;
%         1,0;
%         0,0];
C_yp=[0,1];
C_zp=[0,W_e;
    0,0]; % z=[e;u] = [r-C_yp*x;u]=[-C_yp;0]x + [0;1] u + [1 0;0 0][r;d]
D_zu=[0;W_u];
D_zw=[0,0;
    0,0];
end


%
%  figure;subplot(2,1,1);plot(track(:,1),track(:,[3,2]));
%     set(0,'DefaultTextInterpreter','Latex');set(gca,'FontSize',20);
%     set(gca,'TicklabelInterpreter','Latex');
%     set( get( gca, 'XLabel' ), 'Interpreter', 'Latex' );
%     set( get( gca, 'YLabel' ), 'Interpreter', 'Latex' );
%     title('Tracking performance','FontSize', 24);
%     ylabel('AFR','FontSize', 22);
%    % xlabel('Time [s]','FontSize', 22);
%     L=legend('Output AFR.','Ref. AFR');
%     set(L,'Interpreter','Latex');
%     grid;
%
%     subplot(2,1,2);plot(u(:,1),u(:,2));
%     set(0,'DefaultTextInterpreter','Latex');set(gca,'FontSize',20);
%     set(gca,'TicklabelInterpreter','Latex');
%     set( get( gca, 'XLabel' ), 'Interpreter', 'Latex' );
%     set( get( gca, 'YLabel' ), 'Interpreter', 'Latex' );
%     title('Control effort','FontSize', 24);
%     ylabel('$u(t)$','FontSize', 22);
%     xlabel('Time [s]','FontSize', 22);
%     %L=legend('$J_{th_2}$','$\Vert r_2(t) \Vert_{2,\tau}$');
%     set(L,'Interpreter','Latex');
%     grid;