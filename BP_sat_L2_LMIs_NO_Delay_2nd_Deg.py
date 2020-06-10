import cvxpy as cp
import numpy as np
import mosek


def Sys_Dyn(rho_1, rho_2, delta_1, W_e, W_u): # System Dynamics
    A_p = np.array([
        [-1 / rho_1, 0],
        [-1, -delta_1]])
    B_p = np.array([
        [rho_2/rho_1],
        [0]])
    D_p = np.array([
        [0, .01],
        [1, 0]])  # Omega = [r;d]
    C_y = np.array([
        [1, 0],
        [0, 1]])
    C_z = np.array([
        [0, W_e],
        [0, 0]])  # z = [e;u] = [r - C_yp * x;u]=[-C_yp;0]x + [0;1] u + [1 0;0 0][r;d]
    D_zu = np.array([
        [0],
        [W_u]])
    D_zd = np.array([
        [0, 0],
        [0, 0]])
    return A_p, B_p, D_p, C_y, C_z, D_zu, D_zd


def Mat_No_Delay(rho_1, rho_2, var1, var2, var3, var4, var5): # Matrix polynomial for delay free variables
    output_no_delay = var1 + rho_1 * var2 + 1/2 * (rho_1**2) * var3 + rho_2 * var4 + 1/2 * (rho_2**2) * var5
    return output_no_delay


def Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, var1, var2, var3, var4, var5, var6, var7): # Matrix polynomial for delay-variables
    output_delay = var1 + rho_1 * var2 + rho_2 * var3 + rho_1 * rho_2 * var4 + rho_1_d * var5 + rho_2_d * var6 + rho_1_d*rho_2_d *var7
    return output_delay


# Weights
c_w_1 = .01  # Disturbance weight for the first state
c_w_2 = .01  # Disturbance weight for the second state
C_ew = np.array([
    [c_w_1, 0],
    [0, c_w_2]])
W_e = 1  # Tracking error weight
W_u = 0.01  # Input weight


# Constants
delta_1 = 1e-12  
v_1 = .06

v_2 = 2e-4

theta_max = 8
v_3 = .5
theta_prime = 1

# Disturbance Dynamics
W = np.array([
    [0, .0035*np.pi],
    [-.0035*np.pi, 0]])
V = np.array([
    [1000, 100]])
H = np.array([
    [.01],
    [.01]])
J = 0

# System Dimensions
A_p, B_p, D_p, C_y, C_z, D_zu, D_zd = Sys_Dyn(1, 1, 1, 1, 1)
n_p = A_p.shape[0]
n_d = D_p.shape[1]
n_u = B_p.shape[1]
n_z = C_z.shape[0]
n_y = C_y.shape[0]
n_omega = W.shape[0]
n_d_i = V.shape[0]
n_nu = H.shape[1]

V_bar = np.block([
    [np.zeros((n_u, n_p)), np.zeros((n_u, n_p)), V]])
J_bar = np.block([
    [np.zeros((n_u, n_d)), J]])


# Saturation
u_bar = np.array([
    [120]])
delta = 4e-5  # norm(w(t),2)^-2

# Uncertainty
E = np.array([
    [.02, 0],
    [0, 0]])
G_A = np.array([
    [0.005, 0],
    [0, 0]])
G_B = np.array([
    [0.005],
    [0]])  # zeros(n_p,n_u)
G_D = np.zeros((n_p, n_d))  # zeros(n_p,n_d); np.array([[0.02,0,0],[0,0.02,0.02]])


rho_1_vals = range(100, 201, 50)   # Range of lag time
rho_2_vals = np.array([.45, .55, .66])   # Range of sensitivity
v_1_val = np.array([-v_1, v_1])
v_2_val = np.array([-v_2, v_2])
v_3_val = np.array([-v_3, v_3])

# Tuning Parameters
kappa_vals = np.array([5])  # np.array([.1, 5, 20])
epsilon_vals = np.array([50])  # np.array([1, 50, 200])
ii=0


def main(ii):
    LMI = []
    for epsilon in epsilon_vals:
        for kappa in kappa_vals:
            P_0 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            P_1 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            P_2 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            P_3 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            P_4 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))

            Q_0 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            Q_1 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            Q_2 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            Q_3 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            Q_4 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))

            S_0 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            S_1 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            S_2 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            S_3 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            S_4 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))

            R = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega), PSD=True)
            S1 = cp.Variable((2*n_p + n_omega, 2*n_p + n_omega))
            gamma_sq = cp.Variable(nonneg=True)

            beta = cp.Variable(nonneg=True)

            X_0 = cp.Variable((n_p, n_p))
            X_1 = cp.Variable((n_p, n_p))
            X_2 = cp.Variable((n_p, n_p))
            X_3 = cp.Variable((n_p, n_p))
            X_4 = cp.Variable((n_p, n_p))

            Y_0 = cp.Variable((n_p, n_p))
            Y_1 = cp.Variable((n_p, n_p))
            Y_2 = cp.Variable((n_p, n_p))
            Y_3 = cp.Variable((n_p, n_p))
            Y_4 = cp.Variable((n_p, n_p))

            Z_0 = cp.Variable((n_omega, n_omega))
            Z_1 = cp.Variable((n_omega, n_omega))
            Z_2 = cp.Variable((n_omega, n_omega))
            Z_3 = cp.Variable((n_omega, n_omega))
            Z_4 = cp.Variable((n_omega, n_omega))

            A_K_tilde_0 = cp.Variable((n_p, n_p))
            A_K_tilde_1 = cp.Variable((n_p, n_p))
            A_K_tilde_2 = cp.Variable((n_p, n_p))
            A_K_tilde_3 = cp.Variable((n_p, n_p))
            A_K_tilde_4 = cp.Variable((n_p, n_p))

            B_K_tilde_0 = cp.Variable((n_p, n_y))
            B_K_tilde_1 = cp.Variable((n_p, n_y))
            B_K_tilde_2 = cp.Variable((n_p, n_y))
            B_K_tilde_3 = cp.Variable((n_p, n_y))
            B_K_tilde_4 = cp.Variable((n_p, n_y))

            L_K_tilde_0 = cp.Variable((n_omega, n_p))
            L_K_tilde_1 = cp.Variable((n_omega, n_p))
            L_K_tilde_2 = cp.Variable((n_omega, n_p))
            L_K_tilde_3 = cp.Variable((n_omega, n_p))
            L_K_tilde_4 = cp.Variable((n_omega, n_p))

            L_d_tilde_0 = cp.Variable((n_omega, n_y))
            L_d_tilde_1 = cp.Variable((n_omega, n_y))
            L_d_tilde_2 = cp.Variable((n_omega, n_y))
            L_d_tilde_3 = cp.Variable((n_omega, n_y))
            L_d_tilde_4 = cp.Variable((n_omega, n_y))

            L_y_tilde_0 = cp.Variable((n_omega, n_y))
            L_y_tilde_1 = cp.Variable((n_omega, n_y))
            L_y_tilde_2 = cp.Variable((n_omega, n_y))
            L_y_tilde_3 = cp.Variable((n_omega, n_y))
            L_y_tilde_4 = cp.Variable((n_omega, n_y))

            C_theta_tilde_theta_0 = cp.Variable((n_u, n_p))
            C_theta_tilde_theta_1 = cp.Variable((n_u, n_p))
            C_theta_tilde_theta_2 = cp.Variable((n_u, n_p))
            C_theta_tilde_theta_3 = cp.Variable((n_u, n_p))
            C_theta_tilde_theta_4 = cp.Variable((n_u, n_p))
            C_theta_tilde_theta_5 = cp.Variable((n_u, n_p))
            C_theta_tilde_theta_6 = cp.Variable((n_u, n_p))

            A_theta_tilde_theta_0 = cp.Variable((n_p, n_p))
            A_theta_tilde_theta_1 = cp.Variable((n_p, n_p))
            A_theta_tilde_theta_2 = cp.Variable((n_p, n_p))
            A_theta_tilde_theta_3 = cp.Variable((n_p, n_p))
            A_theta_tilde_theta_4 = cp.Variable((n_p, n_p))
            A_theta_tilde_theta_5 = cp.Variable((n_p, n_p))
            A_theta_tilde_theta_6 = cp.Variable((n_p, n_p))

            B_theta_tilde_theta_0 = cp.Variable((n_p, n_y))
            B_theta_tilde_theta_1 = cp.Variable((n_p, n_y))
            B_theta_tilde_theta_2 = cp.Variable((n_p, n_y))
            B_theta_tilde_theta_3 = cp.Variable((n_p, n_y))
            B_theta_tilde_theta_4 = cp.Variable((n_p, n_y))
            B_theta_tilde_theta_5 = cp.Variable((n_p, n_y))
            B_theta_tilde_theta_6 = cp.Variable((n_p, n_y))

            D_K_theta_0 = cp.Variable((n_u, n_y))
            D_K_theta_1 = cp.Variable((n_u, n_y))
            D_K_theta_2 = cp.Variable((n_u, n_y))
            D_K_theta_3 = cp.Variable((n_u, n_y))
            D_K_theta_4 = cp.Variable((n_u, n_y))
            D_K_theta_5 = cp.Variable((n_u, n_y))
            D_K_theta_6 = cp.Variable((n_u, n_y))

            L_K_theta_tilde_theta_0 = cp.Variable((n_omega, n_p))
            L_K_theta_tilde_theta_1 = cp.Variable((n_omega, n_p))
            L_K_theta_tilde_theta_2 = cp.Variable((n_omega, n_p))
            L_K_theta_tilde_theta_3 = cp.Variable((n_omega, n_p))
            L_K_theta_tilde_theta_4 = cp.Variable((n_omega, n_p))
            L_K_theta_tilde_theta_5 = cp.Variable((n_omega, n_p))
            L_K_theta_tilde_theta_6 = cp.Variable((n_omega, n_p))

            L_y_theta_tilde_theta_0 = cp.Variable((n_omega, n_y))
            L_y_theta_tilde_theta_1 = cp.Variable((n_omega, n_y))
            L_y_theta_tilde_theta_2 = cp.Variable((n_omega, n_y))
            L_y_theta_tilde_theta_3 = cp.Variable((n_omega, n_y))
            L_y_theta_tilde_theta_4 = cp.Variable((n_omega, n_y))
            L_y_theta_tilde_theta_5 = cp.Variable((n_omega, n_y))
            L_y_theta_tilde_theta_6 = cp.Variable((n_omega, n_y))

            E_K_tilde_theta_0 = cp.Variable((n_p, n_u))
            E_K_tilde_theta_1 = cp.Variable((n_p, n_u))
            E_K_tilde_theta_2 = cp.Variable((n_p, n_u))
            E_K_tilde_theta_3 = cp.Variable((n_p, n_u))
            E_K_tilde_theta_4 = cp.Variable((n_p, n_u))
            E_K_tilde_theta_5 = cp.Variable((n_p, n_u))
            E_K_tilde_theta_6 = cp.Variable((n_p, n_u))

            F_K_tilde_theta_0 = cp.Variable((n_omega, n_u))
            F_K_tilde_theta_1 = cp.Variable((n_omega, n_u))
            F_K_tilde_theta_2 = cp.Variable((n_omega, n_u))
            F_K_tilde_theta_3 = cp.Variable((n_omega, n_u))
            F_K_tilde_theta_4 = cp.Variable((n_omega, n_u))
            F_K_tilde_theta_5 = cp.Variable((n_omega, n_u))
            F_K_tilde_theta_6 = cp.Variable((n_omega, n_u))

            T_bar_theta_0 = cp.Variable((n_u, n_u), diag=True)
            T_bar_theta_1 = cp.Variable((n_u, n_u), diag=True)
            T_bar_theta_2 = cp.Variable((n_u, n_u), diag=True)
            T_bar_theta_3 = cp.Variable((n_u, n_u), diag=True)
            T_bar_theta_4 = cp.Variable((n_u, n_u), diag=True)
            T_bar_theta_5 = cp.Variable((n_u, n_u), diag=True)
            T_bar_theta_6 = cp.Variable((n_u, n_u), diag=True)

            G_hat_theta_0 = cp.Variable((n_u, 2*n_p + n_omega))
            G_hat_theta_1 = cp.Variable((n_u, 2*n_p + n_omega))
            G_hat_theta_2 = cp.Variable((n_u, 2*n_p + n_omega))
            G_hat_theta_3 = cp.Variable((n_u, 2*n_p + n_omega))
            G_hat_theta_4 = cp.Variable((n_u, 2*n_p + n_omega))
            G_hat_theta_5 = cp.Variable((n_u, 2*n_p + n_omega))
            G_hat_theta_6 = cp.Variable((n_u, 2*n_p + n_omega))

            G_1_hat_theta_0 = cp.Variable((n_u, 2*n_p + n_omega))
            G_1_hat_theta_1 = cp.Variable((n_u, 2*n_p + n_omega))
            G_1_hat_theta_2 = cp.Variable((n_u, 2*n_p + n_omega))
            G_1_hat_theta_3 = cp.Variable((n_u, 2*n_p + n_omega))
            G_1_hat_theta_4 = cp.Variable((n_u, 2*n_p + n_omega))
            G_1_hat_theta_5 = cp.Variable((n_u, 2*n_p + n_omega))
            G_1_hat_theta_6 = cp.Variable((n_u, 2*n_p + n_omega))

            LMI += [cp.bmat([
                [R, S1],
                [S1.T, R]])>>0]
            LMI += [(delta - beta) >= 0]
            for rho_1 in rho_1_vals:
                for rho_2 in rho_2_vals:
                    A_p, B_p, D_p, C_y, C_z, D_zu, D_zd = Sys_Dyn(rho_1, rho_2, delta_1, W_e, W_u)

                    X = Mat_No_Delay(rho_1, rho_2, X_0, X_1, X_2, X_3, X_4)
                    Y = Mat_No_Delay(rho_1, rho_2, Y_0, Y_1, Y_2, Y_3, Y_4)
                    Z = Mat_No_Delay(rho_1, rho_2, Z_0, Z_1, Z_2, Z_3, Z_4)

                    L_d_tilde = Mat_No_Delay(rho_1, rho_2, L_d_tilde_0, L_d_tilde_1, L_d_tilde_2, L_d_tilde_3, L_d_tilde_4)
                    A_K_tilde = Mat_No_Delay(rho_1, rho_2, A_K_tilde_0, A_K_tilde_1, A_K_tilde_2, A_K_tilde_3, A_K_tilde_4)
                    B_K_tilde = Mat_No_Delay(rho_1, rho_2, B_K_tilde_0, B_K_tilde_1, B_K_tilde_2, B_K_tilde_3, B_K_tilde_4)
                    L_K_tilde = Mat_No_Delay(rho_1, rho_2, L_K_tilde_0, L_K_tilde_1, L_K_tilde_2, L_K_tilde_3, L_K_tilde_4)
                    L_y_tilde = Mat_No_Delay(rho_1, rho_2, L_y_tilde_0, L_y_tilde_1, L_y_tilde_2, L_y_tilde_3, L_y_tilde_4)

                    P_bar = Mat_No_Delay(rho_1, rho_2, P_0, P_1, P_2, P_3, P_4)
                    Q_bar = Mat_No_Delay(rho_1, rho_2, Q_0, Q_1, Q_2, Q_3, Q_4)
                    S_bar = Mat_No_Delay(rho_1, rho_2, S_0, S_1, S_2, S_3, S_4)

                    B_w_hat = cp.bmat([
                        [D_p, (B_p * J)],
                        [(X @ D_p), ((X@B_p)*J)],
                        [(L_d_tilde @ C_y @ D_p), (Z @ H) + ((L_d_tilde @ C_y @ B_p)*J)]])

                    C_hat = cp.bmat([
                        [(C_z @ Y), C_z, (D_zu @ V)],
                        [np.zeros((C_ew.shape[0], n_p)), np.zeros((C_ew.shape[0], n_p)), C_ew]])

                    D_psi = np.block([
                        [-D_zu],
                        [np.zeros((C_ew.shape[0], n_u))]])

                    D_w = np.block([
                        [D_zd, np.dot(D_zu, J)],
                        [np.zeros((C_ew.shape[0], n_d)), np.zeros((C_ew.shape[0], n_nu))]])

                    E_hat = cp.bmat([
                        [E, np.zeros((n_p, n_p)), np.zeros((n_p, n_p))],
                        [(X @ E), np.zeros((n_p, n_p)), np.zeros((n_p, n_p))],
                        [(L_d_tilde @ C_y @ E), np.zeros((n_omega, n_p)), np.zeros((n_omega, n_p))]])

                    for rho_1_d in rho_1_vals:
                        for rho_2_d in rho_2_vals:
                            L_K_theta_tilde_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, L_K_theta_tilde_theta_0,
                                                              L_K_theta_tilde_theta_1, L_K_theta_tilde_theta_2,
                                                              L_K_theta_tilde_theta_3, L_K_theta_tilde_theta_4,
                                                              L_K_theta_tilde_theta_5, L_K_theta_tilde_theta_6)
                            A_theta_tilde_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, A_theta_tilde_theta_0,
                                                            A_theta_tilde_theta_1, A_theta_tilde_theta_2,
                                                            A_theta_tilde_theta_3, A_theta_tilde_theta_4,
                                                            A_theta_tilde_theta_5, A_theta_tilde_theta_6)
                            B_theta_tilde_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, B_theta_tilde_theta_0,
                                                            B_theta_tilde_theta_1, B_theta_tilde_theta_2,
                                                            B_theta_tilde_theta_3, B_theta_tilde_theta_4,
                                                            B_theta_tilde_theta_5, B_theta_tilde_theta_6)
                            C_theta_tilde_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, C_theta_tilde_theta_0,
                                                            C_theta_tilde_theta_1, C_theta_tilde_theta_2,
                                                            C_theta_tilde_theta_3, C_theta_tilde_theta_4,
                                                            C_theta_tilde_theta_5, C_theta_tilde_theta_6)
                            T_bar_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, T_bar_theta_0, T_bar_theta_1,
                                                    T_bar_theta_2, T_bar_theta_3, T_bar_theta_4, T_bar_theta_5,
                                                    T_bar_theta_6)
                            E_K_tilde_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, E_K_tilde_theta_0,
                                                        E_K_tilde_theta_1, E_K_tilde_theta_2, E_K_tilde_theta_3,
                                                        E_K_tilde_theta_4, E_K_tilde_theta_5, E_K_tilde_theta_6)
                            F_K_tilde_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, F_K_tilde_theta_0,
                                                        F_K_tilde_theta_1, F_K_tilde_theta_2, F_K_tilde_theta_3,
                                                        F_K_tilde_theta_4, F_K_tilde_theta_5, F_K_tilde_theta_6)
                            G_hat_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, G_hat_theta_0, G_hat_theta_1,
                                                    G_hat_theta_2, G_hat_theta_3, G_hat_theta_4, G_hat_theta_5,
                                                    G_hat_theta_6)
                            G_1_hat_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, G_1_hat_theta_0, G_1_hat_theta_1,
                                                      G_1_hat_theta_2, G_1_hat_theta_3, G_1_hat_theta_4,
                                                      G_1_hat_theta_5, G_1_hat_theta_6)
                            D_K_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, D_K_theta_0, D_K_theta_1, D_K_theta_2,
                                                  D_K_theta_3, D_K_theta_4, D_K_theta_5, D_K_theta_6)
                            L_y_theta_tilde_theta = Mat_Delay(rho_1, rho_1_d, rho_2, rho_2_d, L_y_theta_tilde_theta_0,
                                                              L_y_theta_tilde_theta_1, L_y_theta_tilde_theta_2,
                                                              L_y_theta_tilde_theta_3, L_y_theta_tilde_theta_4,
                                                              L_y_theta_tilde_theta_5, L_y_theta_tilde_theta_6)
                            Q_bar_delayed = Mat_No_Delay(rho_1_d, rho_2_d, Q_0, Q_1, Q_2, Q_3, Q_4)
                            P_bar_delayed = Mat_No_Delay(rho_1_d, rho_2_d, P_0, P_1, P_2, P_3, P_4)
                            S_h = Mat_No_Delay(rho_1_d, rho_2_d, S_0, S_1, S_2, S_3, S_4)

                            A_d_hat = cp.bmat([
                                [(B_p @ C_theta_tilde_theta), (B_p @ D_K_theta @ C_y), np.zeros((n_p, n_omega))],
                                [A_theta_tilde_theta, (B_theta_tilde_theta @ C_y), np.zeros((n_p, n_omega))],
                                [L_K_theta_tilde_theta, (L_y_theta_tilde_theta @ C_y), np.zeros((n_omega, n_omega))]])
                            K_bb_theta = cp.bmat([
                                [C_theta_tilde_theta, (D_K_theta @ C_y), np.zeros((n_u, n_omega))]])
                            B_psi_hat_theta = cp.bmat([
                                [(-B_p @ T_bar_theta)],
                                [E_K_tilde_theta],
                                [F_K_tilde_theta]])
                            C_d_hat = cp.bmat([
                                [(D_zu @ C_theta_tilde_theta), (D_zu @ D_K_theta @ C_y), np.zeros((n_z, n_omega))],
                                [np.zeros((C_ew.shape[0], n_p)), np.zeros((C_ew.shape[0], n_p)), np.zeros((C_ew.shape[0], n_omega))]])

                            L_13 = R - S1 + A_d_hat
                            L_14 = S1
                            L_15 = B_psi_hat_theta + (G_hat_theta.T) - (K_bb_theta.T) + 2*(V_bar.T)
                            L_16 = B_w_hat
                            L_17 = C_hat.T
                            L_18 = E_hat
                            L_19 = epsilon * cp.bmat([
                                [(G_A @ Y), G_A, (G_B @ V)],
                                [np.zeros((n_p, n_p)), np.zeros((n_p, n_p)), np.zeros((n_p, n_omega))],
                                [np.zeros((n_p, n_p)), np.zeros((n_p, n_p)), np.zeros((n_p, n_omega))]]).T

                            L_22 = -2 * kappa * cp.bmat([
                                [Y, np.eye(n_p), np.zeros((n_p, n_omega))],
                                [np.eye(n_p), X, np.zeros((n_p, n_omega))],
                                [np.zeros((n_omega, n_p)), np.zeros((n_omega, n_p)), Z]]) + (theta_max**2)*R
                            L_23 = kappa*A_d_hat
                            L_24 = np.zeros((2*n_p + n_omega, 2*n_p + n_omega))
                            L_25 = kappa * B_psi_hat_theta
                            L_26 = kappa * B_w_hat
                            L_27 = np.zeros((2*n_p + n_omega, (C_hat.T).shape[1]))
                            L_28 = kappa * E_hat
                            L_29 = np.zeros((2*n_p + n_omega, 3*n_p))

                            L_34 = R - S1.T
                            L_35 = G_1_hat_theta.T + K_bb_theta.T
                            L_36 = np.zeros((2 * n_p + n_omega, B_w_hat.shape[1]))
                            L_37 = C_d_hat.T
                            L_38 = np.zeros((2*n_p + n_omega, 3*n_p))
                            L_39 = epsilon * (cp.bmat([
                                [(G_B @ C_theta_tilde_theta), (G_B @ D_K_theta @ C_y), np.zeros((n_p, n_omega))],
                                [np.zeros((n_p, n_p)), np.zeros((n_p, n_p)), np.zeros((n_p, n_omega))],
                                [np.zeros((n_p, n_p)), np.zeros((n_p, n_p)), np.zeros((n_p, n_omega))]])).T

                            L_44 = -S_h - R
                            L_45 = np.zeros((2*n_p + n_omega, T_bar_theta.shape[1]))
                            L_46 = np.zeros((2*n_p + n_omega, B_w_hat.shape[1]))
                            L_47 = np.zeros((2*n_p + n_omega, (C_hat.T).shape[1]))
                            L_48 = np.zeros((2*n_p + n_omega, 3*n_p))
                            L_49 = np.zeros((2*n_p + n_omega, 3*n_p))

                            L_55 = -4 * T_bar_theta
                            L_56 = 2 * J_bar
                            L_57 = T_bar_theta*(D_psi.T)
                            L_58 = np.zeros((T_bar_theta.shape[0], 3*n_p))
                            L_59 = epsilon * (T_bar_theta @ (np.block([
                                [-G_B],
                                [np.zeros((n_p, n_u))],
                                [np.zeros((n_p, n_u))]])).T)

                            L_66 = -np.eye(B_w_hat.shape[1])
                            L_67 = D_w.T
                            L_68 = np.zeros(((D_w.T).shape[0],3*n_p))
                            L_69 = epsilon * ((np.block([
                                [G_D, np.dot(G_B, J)],
                                [np.zeros((n_p, n_d)), np.zeros((n_p, n_nu))],
                                [np.zeros((n_p, n_d)), np.zeros((n_p, n_nu))]])).T)

                            L_77 = -gamma_sq *np.eye((D_w.T).shape[1])
                            L_78 = np.zeros(((D_w.T).shape[1], 3*n_p))
                            L_79 = np.zeros(((D_w.T).shape[1], 3*n_p))

                            L_88 = -epsilon * np.eye(3*n_p)
                            L_89 = np.zeros((3*n_p, 3*n_p))

                            L_99 = -epsilon * np.eye(3*n_p)

                            for v_1 in v_1_val:
                                for v_2 in v_2_val:
                                    P_dot = v_1 * (P_1 + rho_1 * P_2) + v_2 * (P_3 + rho_2 * P_4)
                                    L_d_tilde_dot = v_1 * (L_d_tilde_1 + rho_1 * L_d_tilde_2) + \
                                                    v_2 * (L_d_tilde_3 + rho_2 * L_d_tilde_4)

                                    A_hat = cp.bmat([
                                        [A_p @ Y, A_p, B_p @ V],
                                        [A_K_tilde, (X @ A_p + B_K_tilde @ C_y), X @ B_p @ V],
                                        [L_K_tilde, (L_d_tilde @ C_y @ A_p) + (L_d_tilde_dot - L_y_tilde) @ C_y, (Z @ W + L_d_tilde @ C_y @ B_p @ V)]])

                                    L_11 = A_hat + A_hat.T + Q_bar + S_bar + P_dot - R
                                    L_12 = P_bar - (cp.bmat([
                                        [Y, np.eye(n_p), np.zeros((n_p, n_omega))],
                                        [np.eye(n_p), X, np.zeros((n_p, n_omega))],
                                        [np.zeros((n_omega, n_p)), np.zeros((n_omega, n_p)), Z]])) + kappa*(A_hat.T)

                                    for v_3 in v_3_val:
                                        L_33 = -2 * R + S1 + S1.T-(1-v_3*theta_prime)*Q_bar_delayed

                                        L = cp.bmat([
                                            [L_11, L_12, L_13, L_14, L_15, L_16, L_17, L_18, L_19],
                                            [L_12.T, L_22, L_23, L_24, L_25, L_26, L_27, L_28, L_29],
                                            [L_13.T, L_23.T, L_33, L_34, L_35, L_36, L_37, L_38, L_39],
                                            [L_14.T, L_24.T, L_34.T, L_44, L_45, L_46, L_47, L_48, L_49],
                                            [L_15.T, L_25.T, L_35.T, L_45.T, L_55, L_56, L_57, L_58, L_59],
                                            [L_16.T, L_26.T, L_36.T, L_46.T, L_56.T, L_66, L_67, L_68, L_69],
                                            [L_17.T, L_27.T, L_37.T, L_47.T, L_57.T, L_67.T, L_77, L_78, L_79],
                                            [L_18.T, L_28.T, L_38.T, L_48.T, L_58.T, L_68.T, L_78.T, L_88, L_89],
                                            [L_19.T, L_29.T, L_39.T, L_49.T, L_59.T, L_69.T, L_79.T, L_89.T, L_99]])
                                        LMI += [L <<0 ]
                                        ii+=1

                            for counter in (range(1, n_u+1)):
                                LMI_sat_d = cp.bmat([
                                    [beta * np.eye(n_u), K_bb_theta[:counter, :] - G_1_hat_theta[:counter, :]],
                                    [(K_bb_theta[:counter, :] - G_1_hat_theta[:counter, :]).T, (u_bar[counter-1]**2)*P_bar_delayed]])
                                LMI += [LMI_sat_d >> 0]

                                LMI_sat = cp.bmat([
                                    [beta * np.eye(n_u), K_bb_theta[:counter, :] - G_hat_theta[:counter, :]],
                                    [(K_bb_theta[:counter, :] - G_hat_theta[:counter, :]).T, (u_bar[counter-1]**2)*P_bar]])
                                LMI += [LMI_sat >> 0]

                    LMI += [Q_bar == Q_bar.T]
                    LMI += [Q_bar >>0]
                    LMI += [S_bar == S_bar.T]
                    LMI += [S_bar >>0]
                    LMI += [P_bar == P_bar.T]
                    LMI += [P_bar_delayed == P_bar_delayed.T]
    prob = cp.Problem(cp.Minimize(gamma_sq), LMI)
    print(ii)
    print(cp.installed_solvers())
    prob.solve(solver=cp.SCS, verbose=True, max_iters=5000, eps=1e-1) # 
    # prob.solve(solver=cp.MOSEK, verbose=True, mosek_params={mosek.dparam.intpnt_co_tol_pfeas: 1.0e-1,
    #                                                         mosek.iparam.intpnt_solve_form: mosek.solveform.primal})
    print("optimal value with:", prob.value)
    print("status:", prob.status)
    print("gamma=", np.sqrt(gamma_sq.value))

main(ii)
