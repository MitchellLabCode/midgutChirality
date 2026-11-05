import numpy as np


## Rodrigues rotation matrix


def rot(vector, angle):
    """3D Rotation matrix from vector, angle"""
    
    c = np.cos(angle)
    s = np.sin(angle)
    K = np.cross(vector, np.identity(vector.shape[0]) * -1)
    
    return np.eye(3) + s*K + (1-c)*K@K


## centerline properties

def _writhe(x):
    """ Compute the writhe of the rod with vertices at positions x (N,3) array 
    Uses algorithm 1a from Klenin & Langowski (2000)
    
    """
    wr = 0
    #print(x.shape)
    N = x.shape[0]
    for i in range(N - 2):
        for j in range(i + 1, N - 1):
            p_1 = x[i, :]
            p_2 = x[i + 1, :]
            p_3 = x[j, :]
            p_4 = x[j + 1, :]

            # Eq 15 in Klenin & Langowski 2000
            r12 = p_2 - p_1
            r34 = p_4 - p_3
            r14 = p_4 - p_1
            r13 = p_3 - p_1
            r23 = p_3 - p_2
            r24 = p_4 - p_2

            #print(r12.shape, r34.shape, r14.shape, r13.shape, r23.shape, r24.shape)

            n1 = np.cross(r13, r14)
            n1 /= np.linalg.norm(n1)
            n2 = np.cross(r14, r24)
            n2 /= np.linalg.norm(n2)
            n3 = np.cross(r24, r23)
            n3 /= np.linalg.norm(n3)
            n4 = np.cross(r23, r13)
            n4 /= np.linalg.norm(n4)

            # Eq 16a in Klenin & Langowski 2000
            omega_star  = (
                np.arcsin(np.dot(n1, n2))
                + np.arcsin(np.dot(n2, n3))
                + np.arcsin(np.dot(n3, n4))
                + np.arcsin(np.dot(n4, n1))
            )

            if np.isnan(omega_star):
                omega_star = 0

            # Eq 16b in Klenin & Langowski 2000
            wr  += omega_star * np.sign(np.dot(np.cross(r34, r12), r13)) / (2 * np.pi)
    return wr


def _curvature(x):
    """Computes the curvature of a rod represented by a Nx3 array of its vertices"""
    kappa = np.zeros(x.shape[0])
    tangent = x[1:,:] - x[:-1,:]
    tangent = tangent / np.linalg.norm(tangent, axis=1)[:,np.newaxis]
    for i in range(1,x.shape[0]-1):
        dots = np.dot(tangent[i,:], tangent[i-1,:])
        if np.abs(np.abs(dots)-1)<1e-7:
            dots = 1.0
        kappa[i] = 2*np.tan(np.arccos(dots)/2.0)
    kappa[0] = kappa[1]
    kappa[-1] = kappa[-2]
    return kappa

def _torsion(x):
    
    tau = np.zeros(x.shape[0]-1) #len N-1
    tangent = x[1:,:] - x[:-1,:] # len N-1
    tangent = tangent / np.linalg.norm(tangent, axis=1)[:,np.newaxis]
    #tau[0] = 0. # because first
    for i in range(1,x.shape[0]-2):
        t1 = tangent[i-1,:]
        t2 = tangent[i,:]
        t3 = tangent[i+1,:]
        axis1 = np.cross(t1, t2)
        axis1 = axis1/ np.linalg.norm(axis1)
        axis2 = np.cross(t2, t3)
        axis2 = axis2/ np.linalg.norm(axis2)
        oAxis = np.cross(t2, axis1)
        cosT = np.dot(axis2, axis1)
        sinT = np.dot(axis2, oAxis)
        tau[i] = np.arctan2(sinT, cosT)
    tau[-1] = tau[-2]
    return tau



def compute_wr(pos):
    """Computes the writhe of the rod
    Input: pos np.array (T, N, 3)
    Output: Wr  np.array (T,)
    """
    if len(pos.shape) == 2:
        return _writhe(pos)
    Wr = np.zeros(pos.shape[0])
    for t in range(pos.shape[0]):
        Wr[t] = _writhe(pos[t,:,:])
    return Wr

def compute_curvature(pos):
    """ compute the curvature
    Input: pos np.array (T, N, 3)
    Output: kappa  np.array (T, N-1)
    """
    if len(pos.shape) == 2:
        return _curvature(pos)
    kappa = np.zeros((pos.shape[0], pos.shape[1]))
    for t in range(pos.shape[0]):
        kappa[t,:] = _curvature(pos[t,:,:])
    return kappa

def compute_torsion(pos):
    """ compute the torsion
    Input: pos np.array (T, N, 3)
    Output: tau  np.array (T, N-2) (?)
    """
    if len(pos.shape) == 2:
        return _torsion(pos)
    tau = np.zeros((pos.shape[0], pos.shape[1]-1))
    for t in range(pos.shape[0]):
        tau[t,:] = _torsion(pos[t,:,:])
    return tau





def compute_coiling(pos):
    """ compute the coiling kappa^2 tau
    Input: pos np.array (T, N, 3)
    Output: khi  np.array (T,N-2)
    """
    if len(pos.shape) == 2:
        return _curvature(pos)
    khi = np.zeros((pos.shape[0], pos.shape[1]))
    for t in range(pos.shape[0]):
        khi[t,:] = _curvature(pos[t,:,:]) * _torsion(pos[t,:,:])
    return khi


def compute_intcoil(pos):
    
    return

def compute_perversion(pos):
    """ derivative of the coiling kappa^2 tau """
    
    return

def _sphereAreaTri(a, b, c):
    """ Computes the spherical area spanned by the spherical triangle (a,b,c), where a, b,c are 3d np arrays."""
    return 2* np.arctan( np.dot(np.cross(a,b), c) / (1+ np.dot(a,b) + np.dot(b,c) + np.dot(a,c)))

def _sAreaPos(pos):
    """ Computes the spherical area spanned by the sequence of points in the Nx3 array of vectors pos.
    vectors are normalized in pre-processing"""
    
    if pos.shape[0] < 3:
        raise Exception("Input too short")
    pos_norm = pos / np.norm(pos, axis=2)[:,np.newaxis]
    s = 0.
    for k in range(1,pos.shape[0]-1):
        s += _sphereAreaTri(pos_norm[0,:], pos_norm[k,:], pos_norm[k+1,:])
    return s
    

def quasi_local(pos, w):
    """ Quasi-local coiling measure with width 2w"""
    q = np.zeros(pos.shape[0] - 2)
    for k in range(0, w):
        q[k] = _sAreaPos(pos[0:k+w,:])
    for k in range(pos.shape[0] - w, pos.shape[0] - 2):
        q[k] = _sAreaPos(pos[k-w:-1,:])
    for k in range(w, pos.shape[0] - w):
        q[k] = _sAreaPos(pos[k-w:k+w,:])
    return q
    

## Framed curve properties


def compute_twist(center_line, d1):
    """
    Compute the twist of a rod using center_line and d1.
    Methods used in this function are adapted from method 2a Klenin & Langowski 2000 paper.

    Parameters
    ----------
    center_line : ndarray (N, 3)
                rod node positions.
    d1 : ndarray (N-1, 3)
                        rod elements normal direction.

    Returns
    -------
    total_twist : float
    local_twist : ndarray (N-2,)
    """
    N = center_line.shape[0]

    total_twist = 0.0
    local_twist = np.zeros(N - 2)

    # Compute s vector
    # s is the tangent

    #print("center_line.shape, d1.shape:", center_line.shape, d1.shape)

    s = center_line[1:] - center_line[:-1]

    # Compute tangents unit directors
    tangent = s / np.linalg.norm(s, axis=1, keepdims=True)

    #print('tangent shape:', tangent.shape)
    #print('d1[:-1,:] shape:', d1[:-1,:].shape)
    if d1.shape[0] == N - 1:
        d1_ = d1
    else:
        d1_ = d1[:-1,:]

    # Compute the projection of normal collection (d1) on normal-binormal plane
    projection_of_normal_collection = d1_ - \
        (tangent * d1_).sum(axis=1, keepdims=True) * tangent
    projection_of_normal_collection = projection_of_normal_collection / \
        np.linalg.norm(projection_of_normal_collection, axis=1, keepdims=True)

    # Eq 27 in Klenin & Langowski 2000
    # p is defined on interior nodes
    p = np.cross(s[:-1], s[1:])
    p = p / np.linalg.norm(p, axis=1, keepdims=True)

    # Compute the angle we need to turn d1 around s to get p
    # Sign part tells whether d1 must be rotated ccw(+) or cw(-) around s
    alpha = np.sign((np.cross(projection_of_normal_collection[:-1], p) * s[:-1]).sum(axis=1)) * \
        np.arccos((projection_of_normal_collection[:-1] * p).sum(axis=1))

    gamma = np.sign((np.cross(p, projection_of_normal_collection[1:]) * s[1:]).sum(axis=1)) * \
        np.arccos((projection_of_normal_collection[1:] * p).sum(axis=1))

    # Normalize angles to a full rotation
    alpha = alpha / (2 * np.pi)
    gamma = gamma / (2 * np.pi)

    twist_temp = alpha + gamma

    # Ensure twist is between (-1/2 to 1/2)
    twist_temp[twist_temp > 0.5] -= 1
    twist_temp[twist_temp < -0.5] += 1

    # Check for NaN values (occur when tangents are parallel)
    twist_temp[np.isnan(twist_temp)] = 0

    # Store twist values
    local_twist[:] = twist_temp
    total_twist = np.sum(twist_temp)

    return total_twist, local_twist

def excess_twist():
    
    return
    
        
def plot_tube(x,Q, color):
    
    return

def skewMat2Vec(M):
    return np.array([M[2,1], - M[2,0], M[1,0]])

def skewVec2Mat(u):
    return np.array([[0., -u[2], u[1]],
                     [u[2], 0., -u[0]],
                     [-u[1], u[0], 0.]])

def logrot(M):
    """ matrix logarithm """
    u = 0.5 * skewMat2Vec(M - M.T)
    if np.linalg.norm(u) == 0.:
        return np.zeros(3)
    else:
        u /= np.linalg.norm(u)
        theta = np.arctan2(-np.trace(skewVec2Mat(u) @  M )/2, 0.5*(np.trace(M) - 1.0))
    return theta*u

def compute_omega(Qs_arr):
    """
    Compute the rotation vector from the frames at different time points
    parameters:
        Qs_arr : (num_times x num_space x 3 x 3) array
    returns:
        omega (num_times -1) x num_space x 3 array
    """
    num_t, num_s, _, _ = Qs_arr.shape
    omega = np.zeros((num_t-1, num_s, 3))
    for t in range(num_t-1):
        for i in range(num_s):
            omega[t,i,:] = logrot(np.squeeze(Qs_arr[t,i,:,:]).T @ np.squeeze(Qs_arr[t+1,i,:,:]))
    return omega

def compute_voronoi(x):
    
    segments_l = np.linalg.norm(x[:,1:,:] - x[:,:-1,:],axis=2)
    return 0.5*(segments_l[:,1:] + segments_l[:,:-1])

def compute_darboux(Qs_arr, x):
    """
    Compute the Darboux vector from the frames at different time points
    parameters:
        Qs_arr : (num_times x num_space x 3 x 3) array
        x : (num_times x (num_space+1) x 3 x 3) array
    returns:
        kappa (num_times x (num_space-1) x 3 array
    """
    num_t, num_s, _, _ = Qs_arr.shape
    kappa = np.zeros((num_t, num_s-1, 3))
    voronois = compute_voronoi(x)
    for t in range(num_t):
        for i in range(num_s-1):
            kappa[t,i,:] = logrot(np.squeeze(Qs_arr[t,i,:,:]).T @ np.squeeze(Qs_arr[t,i+1,:,:])) / voronois[t,i]
    return kappa
    