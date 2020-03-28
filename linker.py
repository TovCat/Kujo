import reader_cif as rc
import reader_orca as ro
import numpy as np

def link(cluster: rc.Cluster, orca_mol: rc.Molecule, mu: np.array(1, 3)):

    def q_mul(q1: np.array(1, 4), q2: np.array(1, 4)):
        mul = np.zeros(1, 4)
        mul[0, 0] = q1[0, 0]*q2[0, 0] - q1[0, 1]*q2[0, 1] - q1[0, 2]*q2[0, 2] - q1[0, 3]*q2[0, 3]
        mul[0, 1] = q1[0, 0]*q1[0, 1] + q1[0, 1]*q2[0, 0] + q1[0, 2]*q2[0, 3] - q1[0, 3]*q2[0, 2]
        mul[0, 2] = q1[0, 0]*q2[0, 2] - q1[0, 1]*q2[0, 3] + q1[0, 2]*q2[0, 0] + q1[0, 3]*q2[0, 1]
        mul[0, 3] = q1[0, 0]*q2[0, 3] + q1[0, 1]*q2[0, 2] - q1[0, 2]*q2[0, 1] + q1[0, 3]*q2[0, 0]
        return mul

    mu_list = []
    for i in range(len(cluster.molecules)):
        cluster.molecules[i].inertia()
        norm = np.cross(cluster.molecules[i].inertia_eig_vec[0:1], orca_mol.inertia_eig_vec[0:1])
        theta_1 = np.arccos(np.dot(cluster.molecules[i].inertia_eig_vec[0:1], orca_mol.inertia_eig_vec[0:1]))
        mu_q = np.array([0, mu[0, 0], mu[0, 1], mu[0, 2]])
        q1 = np.array([np.cos(theta_1/2), norm[0, 0] * np.sin(theta_1/2), norm[0, 1] * np.sin(theta_1/2),
                       norm[0, 2] * np.sin(theta_1/2)])
        y_q = np.array([0, cluster.molecules[i].inertia_eig_vec[0, 0], cluster.molecules[i].inertia_eig_vec[0, 1],
                        cluster.molecules[i].inertia_eig_vec[0, 2]])
        y_rot = q_mul(q_mul(q1, y_q), q1)
        theta_2 = np.arccos(np.dot(cluster.molecules[i].inertia_eig_vec[1:2], y_rot))
        q2 = np.array([np.cos(theta_2/2), cluster.molecules[i].inertia_eig_vec[0, 0] * np.sin(theta_2/2),
                       cluster.molecules[i].inertia_eig_vec[0, 1] * np.sin(theta_2/2),
                       cluster.molecules[i].inertia_eig_vec[0, 2] * np.sin(theta_2/2)])
        q_final = q_mul(q1, q2)
        mu_q = q_mul(q_mul(q_final, mu_q), q_final)
        mu_rotated = np.array([mu_q[0, 1], mu_q[0, 2], mu_q[0, 3]])
        mu_list.append(mu_rotated)