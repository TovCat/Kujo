import reader_cif as rc
import reader_orca as ro
import numpy as np

def link(cluster: rc.Cluster, orca_mol: rc.Molecule, mu: np.array):

    def q_mul(quaternion0, quaternion1):
        w0, x0, y0, z0 = quaternion0
        w1, x1, y1, z1 = quaternion1
        return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                         x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                         -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                         x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)

    mu_list = []
    for i in range(len(cluster.molecules)):
        cluster.molecules[i].inertia()
        norm = np.cross(cluster.molecules[i].inertia_eig_vec[0:1], orca_mol.inertia_eig_vec[0:1])
        theta_1 = np.arccos(np.dot(np.squeeze(cluster.molecules[i].inertia_eig_vec[0:1]),
                                   np.squeeze(orca_mol.inertia_eig_vec[0:1])))
        mu_q = np.array([0, mu[0, 0], mu[0, 1], mu[0, 2]])
        q1 = np.array([np.cos(theta_1/2), norm[0, 0] * np.sin(theta_1/2), norm[0, 1] * np.sin(theta_1/2),
                       norm[0, 2] * np.sin(theta_1/2)])
        y_q = np.array([0, cluster.molecules[i].inertia_eig_vec[0, 0], cluster.molecules[i].inertia_eig_vec[0, 1],
                        cluster.molecules[i].inertia_eig_vec[0, 2]])
        q_temp = q_mul(q1, y_q)
        y_rot_q = q_mul(q_temp, q1)
        y_rot = np.array([y_rot_q[1], y_rot_q[2], y_rot_q[3]])
        theta_2 = np.arccos(np.dot(np.squeeze(cluster.molecules[i].inertia_eig_vec[1:2]), np.squeeze(y_rot)))
        q2 = np.array([np.cos(theta_2/2), cluster.molecules[i].inertia_eig_vec[0, 0] * np.sin(theta_2/2),
                       cluster.molecules[i].inertia_eig_vec[0, 1] * np.sin(theta_2/2),
                       cluster.molecules[i].inertia_eig_vec[0, 2] * np.sin(theta_2/2)])
        q_final = q_mul(q1, q2)
        q_temp = q_mul(q_final, mu_q)
        mu_q = q_mul(q_temp, q_final)
        mu_rotated = np.array([mu_q[1], mu_q[2], mu_q[3]])
        mu_list.append(mu_rotated)

    return mu_list