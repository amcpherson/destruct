import itertools
import numpy as np
import scipy.stats as st

# High copy probability
p_hc = 1e-16

# Maximum copy deviation
max_dev = 2.

# Deviation probability
p_dev = 0.1

def create_local_grid(a, b):
    mat, pat = np.meshgrid(np.arange(max(0.0, a - 2.0), max(1.0, a + 3.0)), np.arange(max(0.0, b - 2.0), max(1.0, b + 3.0)))
    return zip(mat.flatten(), pat.flatten())

mat_pat_cached = dict()
def create_mat_pat(x, y, s, t):
    a, b = np.rint((x - s) / t), np.rint((y - s) / t)
    if (a, b) in mat_pat_cached:
        return mat_pat_cached[(a, b)]
    mat_pat = list(set(create_local_grid(a, b) + create_local_grid(b, a)))
    mat = np.array([mp[0] for mp in mat_pat])
    pat = np.array([mp[1] for mp in mat_pat])
    mat_pat_cached[(a, b)] = (mat, pat)
    return mat, pat

def resp(x, y, s, t, v_x, v_y, mat, pat):
    try:
        w_0 = -np.square(x - (s + mat * t)) / (2.0 * v_x) - np.square(y - (s + pat * t)) / (2.0 * v_y)
        w_1 = -np.square(x - (s + pat * t)) / (2.0 * v_x) - np.square(y - (s + mat * t)) / (2.0 * v_y)
        w_0 = np.exp(w_0 - w_0.max()) / np.sum(np.exp(w_0 - w_0.max()))
        w_1 = np.exp(w_1 - w_1.max()) / np.sum(np.exp(w_1 - w_1.max()))
        return np.array([w_0, w_1])
    except ValueError:
        print (x, y, s, t, v_x, v_y, mat, pat)
        raise

def calc_n1(x, y, r, mat, pat):
    return np.sum((r[0] + r[1]))

def calc_n2(x, y, r, mat, pat):
    return (x + y) * np.sum((r[0] + r[1]))

def calc_n3(x, y, r, mat, pat):
    return np.sum((r[0] + r[1]) * (mat + pat))

def calc_n4(x, y, r, mat, pat):
    return np.sum((r[0] + r[1]) * (np.square(mat) + np.square(pat)))

def calc_n5(x, y, r, mat, pat):
    return np.sum((r[0] * (x * mat + y * pat) + r[1] * (x * pat + y * mat)))

def optimize_s_t(X, Y, L, s, t, v_x, v_y, num_iter=20):
    for _ in xrange(num_iter):
        n1 = 0.0
        n2 = 0.0
        n3 = 0.0
        n4 = 0.0
        n5 = 0.0
        for x, y, l in zip(X, Y, L):
            mat, pat = create_mat_pat(x, y, s, t)
            r = resp(x, y, s, t, v_x, v_y, mat, pat)
            n1 += l * calc_n1(x, y, r, mat, pat)
            n2 += l * calc_n2(x, y, r, mat, pat)
            n3 += l * calc_n3(x, y, r, mat, pat)
            n4 += l * calc_n4(x, y, r, mat, pat)
            n5 += l * calc_n5(x, y, r, mat, pat)
        s, t = np.linalg.solve(np.array([[2*n1,n3], [n3,n4]]), np.array([n2,n5]))
    return s, t

def xy_prob_z0(x, y, s, t, v_x, v_y):
    mat, pat = create_mat_pat(x, y, s, t)
    w_0 = -np.square(x - (s + mat * t)) / (2.0 * v_x) - np.square(y - (s + pat * t)) / (2.0 * v_y)
    w_1 = -np.square(x - (s + pat * t)) / (2.0 * v_x) - np.square(y - (s + mat * t)) / (2.0 * v_y)
    prob = np.sum(np.exp(w_0)) + np.sum(np.exp(w_1))
    return prob

def xy_prob_z1(x, y, s, t, f, v_x, v_y):
    mat, pat = create_mat_pat(x, y, s, t)
    prob = 0.0
    for g_m in (0, 1):
        for g_p in (0, 1):
            for d_m in np.arange(0., max_dev+1.):
                for d_p in np.arange(0., max_dev+1.):
                    if d_m == 0. and d_p == 0.:
                        continue
                    fg_m = g_m * f + (1-g_m) * (1-f)
                    fg_p = g_p * f + (1-g_p) * (1-f)
                    mat_d = (1-fg_m)*mat + fg_m*(mat+d_m)
                    pat_d = (1-fg_p)*pat + fg_p*(pat+d_p)
                    p_d_mp = d_m * np.log(p_dev) + d_p * np.log(p_dev)
                    w_0 = -np.square(x - (s + mat_d * t)) / (2.0 * v_x) - np.square(y - (s + pat_d * t)) / (2.0 * v_y) + p_d_mp
                    w_1 = -np.square(x - (s + pat_d * t)) / (2.0 * v_x) - np.square(y - (s + mat_d * t)) / (2.0 * v_y) + p_d_mp
                    prob += np.sum(np.exp(w_0))
                    prob += np.sum(np.exp(w_1))
    return prob

def z_posterior(x, y, s, t, f, v_x, v_y):
    z_post = np.array([xy_prob_z0(x, y, s, t, v_x, v_y), xy_prob_z1(x, y, s, t, f, v_x, v_y)])
    return z_post / np.sum(z_post)

def z_sample(x, y, s, t, f, v_x, v_y):
    z_post = z_posterior(x, y, s, t, f, v_x, v_y)
    try:
        return np.random.choice(len(z_post), p=z_post)
    except ValueError:
        print x, y, s, t, f, v_x, v_y, z_post
        raise

def max_major_minor_posterior(x, y, s, t, f, v_x, v_y):
    mat, pat = create_mat_pat(x, y, s, t)
    best_prob = list()
    best_copies = list()
    for g_m, g_p, d_m, d_p in itertools.product((0, 1), (0, 1), np.arange(0., max_dev+1.), np.arange(0., max_dev+1.)):
        fg_m = g_m * f + (1-g_m) * (1-f)
        fg_p = g_p * f + (1-g_p) * (1-f)
        mat_d = (1-fg_m)*mat + fg_m*(mat+d_m)
        pat_d = (1-fg_p)*pat + fg_p*(pat+d_p)
        p_d_mp = d_m * np.log(p_dev) + d_p * np.log(p_dev)
        w_0 = -np.square(x - (s + mat_d * t)) / (2.0 * v_x) - np.square(y - (s + pat_d * t)) / (2.0 * v_y) + p_d_mp
        w_1 = -np.square(x - (s + pat_d * t)) / (2.0 * v_x) - np.square(y - (s + mat_d * t)) / (2.0 * v_y) + p_d_mp
        prob = np.exp(w_0) + np.exp(w_1)
        best = np.argmax(prob)
        best_prob.append(prob[best])
        b_m = (mat[best]+d_m, mat[best])[g_m]
        b_p = (pat[best]+d_p, pat[best])[g_p]
        d_m = (-d_m, d_m)[g_m]
        d_p = (-d_p, d_p)[g_p]
        if b_m > b_p:
            best_copies.append((b_m, b_p, d_m, d_p))
        else:
            best_copies.append((b_p, b_m, d_p, d_m))
    best = np.argmax(np.array(best_prob))
    major1 = (x - s) / t
    minor1 = (y - s) / t
    if abs(best_copies[best][0] - np.rint(major1)) > 1 or abs(best_copies[best][1] - np.rint(minor1)) > 1:
        print 'warning: abnormal region'
        print best_copies[best], major1, minor1
        print (x, y, s, t, f, v_x, v_y)
    return best_copies[best]

def estimate_model_prob(major, minor, length, s, t, v_major, v_minor):

    f = 0.25

    model_probs = list()
    ss = list()
    ts = list()
    fs = list()
    avg_z = np.zeros(len(major))

    burnin = 5
    num_iter = 100

    z = np.zeros(len(major))
    
    print s, t

    for i in xrange(burnin + num_iter):

        # Optimize s and t based on dominant regions
        s, t = optimize_s_t(major[z==0], minor[z==0], length[z==0], s, t, v_major, v_minor, num_iter=5)

        # Sample assignment to clusters from posterior over clusters
        z_new = np.zeros(len(major))
        for idx in np.random.permutation(len(major)):
            z_new[idx] = z_sample(major[idx], minor[idx], s, t, f, v_major, v_minor)
        z = z_new

        # Max likelihood estimate of f
        major_frac = np.absolute(np.rint((major[z==1] - s) / t) - ((major[z==1] - s) / t))
        minor_frac = np.absolute(np.rint((minor[z==1] - s) / t) - ((minor[z==1] - s) / t))
        frac = np.maximum(major_frac, minor_frac)
        f = frac.mean()
        
        # Revert back to initial value for empty cluster
        if np.isnan(f):
            f = 0.25

        # Calculate model probability given current s, t and z
        model_prob = 0.0
        for idx in xrange(len(major)):
            if z[idx] == 0:
                p_x = xy_prob_z0(major[idx], minor[idx], s, t, v_major, v_minor)
            else:
                p_x = xy_prob_z1(major[idx], minor[idx], s, t, f, v_major, v_minor)
            model_prob += length[idx] * np.log(max(p_x, p_hc))

        print s, t, np.sum(z), len(major), model_prob

        model_probs.append(model_prob)
        ss.append(s)
        ts.append(t)
        fs.append(f)

        if i >= burnin:
            avg_z += z

    model_probs = np.array(model_probs)
    ss = np.array(ss)
    ts = np.array(ts)
    fs = np.array(fs)

    avg_z /= float(num_iter)

    best_major_minor = list()
    for idx in xrange(len(major)):
        best_major_minor.append(max_major_minor_posterior(major[idx], minor[idx], s, t, fs[model_probs.argmax()], v_major, v_minor))
    best_major_minor = np.array(best_major_minor).T

    return model_probs.max(), ss[model_probs.argmax()], ts[model_probs.argmax()], fs[model_probs.argmax()], avg_z, best_major_minor

