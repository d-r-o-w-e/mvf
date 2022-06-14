import he_mesh as hm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc
import random
import dionysus as d
from itertools import chain, combinations
import os
import subprocess
import glob

# from the itertools page
def powerset(iterable):
    "powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

def D(M, v, mu, eps):
    # takes in a vector field (function) v on M, and returns v_bar: V -> M
    # v_bar is a "snapping" of the vectors onto the simplices of M
    # mu defines an angle threshold, and eps defines a length threshold
    
    v_bar = {}
    
    for vtx in M.V:
        
        # first, evaluate v at the vertex
        v_x = v(vtx.pos)
        
        adj_hes = []
        adj_fs = []
        
        # get the vertex normal (face normals weighted by area).  assumes each vertex is adj to only 1 boundary
        n = np.zeros(3)
        he0 = vtx.he
        he = he0
        i = 0
        while not he is he0 or i == 0:
            # collect all the adjacent hes for later
            adj_hes += [he]
            adj_fs += [he.face]
            
            # don't sum for boundary elements
            if not he.face is None:
                v0 = he.next.vert.pos - he.vert.pos
                v1 = he.next.next.vert.pos - he.vert.pos
                n += 0.5*np.cross(v0, v1)

                # go to next element of star (?)
                he = he.next.next.twin
            else:
                # if we're on the boundary, loop backward until we aren't
                he1 = he
                while 1:
                    if he1.twin.face is None:
                        he = he1
                        break
                    else:
                        he1 = he1.twin.next
            i += 1
        
        # normalise n
        n = n / np.linalg.norm(n)
        
        # now, get any edge vec e0:
        e0 = vtx.he.next.vert.pos - vtx.pos
        # generate an orthonormal basis with the corrected e0 (basis is e0, e1, n)
        b0 = e0 - np.dot(e0, n)*n
        b0 = b0 / np.linalg.norm(b0)
        b1 = np.cross(n, b0)
        
        # assert that all are length 1
        assert(abs(np.linalg.norm(b0) - 1) < 1e-6)
        assert(abs(np.linalg.norm(b1) - 1) < 1e-6)
        assert(abs(np.linalg.norm(n) - 1) < 1e-6)
        
        # get vectors representing each of the hes and faces adjacent to the vertex
        adj_hes_vec = np.array([h.next.vert.pos - h.vert.pos for h in adj_hes]) # k x 3
        adj_fs_vec = np.array([f.center() - vtx.pos if not f is None else np.zeros(3) for f in adj_fs]) # k x 3
        
        # get their dot products with b0 and b1:
        basis = np.vstack([b0, b1]) # 2 x 3
        he_reps = basis @ adj_hes_vec.T # 2 x k
        f_reps = basis @ adj_fs_vec.T # 2 x k
        v0 = basis @ v_x[:, None] # 2 x 1.  name collision sorry, shouldn't cause issues
        
        if np.linalg.norm(v0) < eps:
            v_bar[vtx] = vtx # point to self
            continue
        
        # normalise all of them
        he_reps = he_reps / np.linalg.norm(he_reps, axis=0)
        f_reps = f_reps / np.where(np.abs(np.linalg.norm(f_reps, axis=0)) < 1e-6, np.ones_like(f_reps), np.linalg.norm(f_reps, axis=0))
        v0 = v0 / np.linalg.norm(v0, axis=0)
        
        assert(np.all(np.linalg.norm(he_reps, axis=0) - 1 < 1e-6))
        assert(np.all(np.linalg.norm(f_reps, axis=0) - 1 < 1e-6))
        assert(np.all(np.linalg.norm(v0, axis=0) - 1 < 1e-6))
        
        # get dot products with all those vectors
        a_he = np.argmax((v0.T @ he_reps)[0, :])
        m_he = np.clip(np.amax((v0.T @ he_reps)[0, :]), -1, 1)
        a_f = np.argmax((v0.T @ f_reps)[0, :])
        
        # determine whether the vector is close enough to an edge to snap to it
        # also snap to nearest edge if face is None
        if np.arccos(m_he) < mu or adj_fs[a_f] is None:
            v_bar[vtx] = adj_hes[a_he].edge
            continue
        
        # otherwise, snap to nearest face
        v_bar[vtx] = adj_fs[a_f]
        
    return v_bar

def CMVF(M, v, mu, eps, outer_ring=False):
    # takes in a mesh and a vector field to be evaluated at the vertices.  
    # returns theta, a function whose preimages are the multivectors
    # if outer_ring is true, make a cycle on the outside
    
    # discretise v.  v_bar: V -> M
    v_bar = D(M, v, mu, eps)
    
    # put a cycle around the outside if outer_ring
    if outer_ring:
        for h in m.HE:
            if h.face is None:
                v_bar[h.vert] = h.edge
    
    # initialise theta
    theta = {}
    
    # SINGLETONS
    # for all simplices, set their theta to themselves
    for v in M.V:
        theta[v] = v
    for e in M.E:
        theta[e] = e
    for f in M.F:
        theta[f] = f
    
    # EDGES
    # if both endpoints of an edge point the same way, pair that edge with the face
    for e in M.E:
        x_m = e.he.vert
        x_p = e.he.twin.vert
        
        # get x_m and x_p as vectors
        vec_m = v_bar[x_m].center() - x_m.pos
        vec_p = v_bar[x_p].center() - x_p.pos
        
        # get an axis perpendicular to e
        face_adj = e.he.face
        if face_adj is None:
            face_adj = e.he.twin.face
        e_perp = face_adj.center() - e.center()
        
        # normalise them all
        vec_m = vec_m / np.where(np.linalg.norm(vec_m) == 0, 1, np.linalg.norm(vec_m))
        vec_p = vec_p / np.where(np.linalg.norm(vec_p) == 0, 1, np.linalg.norm(vec_p))
        e_perp = e_perp / np.linalg.norm(e_perp)
        
        # project both endpoints onto e_perp; if greater than 0, pair edge with that face
        proj_m = np.dot(vec_m, e_perp)
        proj_p = np.dot(vec_p, e_perp)
        if proj_m > 0 and proj_p > 0:
            theta[e] = face_adj
        elif proj_m < 0 and proj_m < 0 and not e.he.twin.face is None:
            theta[e] = e.he.twin.face
        
    # POINTS
    for v in M.V:
        # if v_bar pointing to self, continue
        if v_bar[v] is v: 
            continue
        
        # if paired with edge, connect to that edge. 
        elif isinstance(v_bar[v], hm.Edge):
            theta[v] = theta[v_bar[v]]
        
        # if paired with a face:
        elif isinstance(v_bar[v], hm.Face):
            # get the two edges adjacent to it on that face, e1 and e2
            f = v_bar[v]
            he1 = f.he
            # cycle until we get the right halfedges.  runs at most 3 times
            while he1.vert != v:
                he1 = he1.next
            he2 = he1.next.next # the other halfedge on that face
            
            e1 = he1.edge
            e2 = he2.edge
            
            # determine whether or not e1 and e2 are singletons right now
            t1 = theta[e1] == e1
            t2 = theta[e2] == e2
            
            # get the other endpoints of e1 and e2
            x1s1 = he1.next.vert
            x2s2 = he2.vert
            
            # get e1 (e2) or one of its adjacent faces
            opp1 = filter(lambda x: not x is None, [e1, e1.he.face, e1.he.twin.face])
            opp2 = filter(lambda x: not x is None, [e2, e2.he.face, e2.he.twin.face])
            
            # determine whether or not to pair the vertex with a face / edge
            c1 = t1 and (v_bar[x1s1] in opp1)
            c2 = t2 and (v_bar[x2s2] in opp2)
            
            if c1 and c2:
                continue
            if not c1 and not c2:
                theta[e1] = f
                theta[e2] = f
                theta[v] = f
            else:
                if c1 and t2:
                    theta[v] = e1
                if c2 and t1:
                    theta[v] = e2
    return theta

def cl(A):
    # A is a set of simplices.  cl(A), the closure of A, is all faces of the simplices in A
    clA = A.copy()
    # add every edge of every face
    to_add = set()
    for f in clA:
        if isinstance(f, hm.Face): 
            to_add |= {f.he.edge, f.he.next.edge, f.he.next.next.edge}
    clA |= to_add
    # add every vertex of every edge
    to_add = set()
    for e in clA:
        if isinstance(e, hm.Edge):
            to_add |= {e.he.vert, e.he.twin.vert}
    clA |= to_add
    # that's it.  return
    return clA
        
def mo(A):
    # A is a set of simplices.  mo(A), the mouth of A, is cl(A)\A
    return cl(A) - A

def pf(S, theta):
    # given a set of simplices S, push it forward (explore everywhere reachable) using DFS
    # theta: M -> M takes elements in the same multivector to the same element.  the output of CMVF
    
    # preim: M -> P(M) takes theta's images to its preimages (elements in the multivector)
    # preim[theta[sigma]] will take an element to a list of its multivector elements
    # get preimage of theta
    preim = {}
    for x in theta.keys():
        if theta[x] in preim.keys():
            preim[theta[x]].add(x)
        else:
            preim[theta[x]] = {x}
    
    pfS = S.copy() # S is trivially a part of its own pf
    
    # iterative DFS (since the meshes may be huge)
    # mirrors pseudocode from https://en.wikipedia.org/wiki/Depth-first_search
    visited = set()
    stack = [s for s in S]
    while len(stack) > 0:
        sigma = stack[-1] # read the last element
        stack = stack[:-1] # pop the last element
        if not sigma in visited:
            # visit sigma
            visited.add(sigma)
            pfS.add(sigma)
            # get all the neighbors of sigma wrt the dynamics 
            # do not add sigma itself, it's already in pfS!
            
            # the multivector of sigma \ sigma
            mv = preim[theta[sigma]] - {sigma}
            # the closure of sigma, minus sigma (mouth)
            clsigma = cl({sigma}) - {sigma}
            # neighbors is [sigma]_V union cl(sigma)
            neighbors = mv | clsigma
            # push all the neighbors to the stack and continue
            stack += list(neighbors)
            
    return pfS

def is_critical(sigma):
    # returns True if the mv is critical, False otherwise
    # critical := H_p (cl(A), mo(A)) != 0 in some dimension p
    # assume sigma is a set of vertices, edges, faces
    
    closure = [d.Simplex(a.simplex()) for a in cl(sigma)]
    mouth = [d.Simplex(a.simplex()) for a in mo(sigma)]
    
    f = d.Filtration(closure)
    f0 = d.Filtration(mouth)
    m = d.homology_persistence(f, relative=f0)
    dgms = d.init_diagrams(m, f)
    for dgm in dgms:
        if dgm:
            return True
    return False
            
def is_regular(sigma):
    return not is_critical(sigma)

def inv(theta):
    # compute the invariant part of the mesh.  
    # include all critical multivectors, as well as any multivector which can return to itself in the mv field
    invN = set()
    
    # get the preimage
    preim = {}
    for x in theta.keys():
        if theta[x] in preim.keys():
            preim[theta[x]].add(x)
        else:
            preim[theta[x]] = {x}
            
    # all critical multivectors are part of invN
    for x in preim.keys():
        # print(x)
        if is_critical(preim[x]):
            invN |= preim[x]
        else:
            for x in preim.keys():
                mv = preim[x]
                
                visited = set()
                stack = [s for s in mo(mv)]
                while len(stack) > 0:
                    sigma = stack[-1] # read the last element
                    stack = stack[:-1] # pop the last element
                    if not sigma in visited:
                        
                        # if we've returned to the same multivector, add that multivector to invN
                        if sigma in mv:
                            invN |= mv
                            break
                        
                        # visit sigma
                        visited.add(sigma)
                        # get all the neighbors of sigma wrt the dynamics 
                        # do not add sigma itself!
                        curr_mv = preim[theta[sigma]] - {sigma}
                        # the closure of sigma, minus sigma (mouth)
                        mosigma = cl({sigma}) - {sigma}
                        # neighbors is [sigma]_V union cl(sigma)
                        neighbors = curr_mv | mosigma
                        # push all the neighbors to the stack and continue
                        stack += list(neighbors)
                        
    return invN

def conley_index_persistence(v, m, times):
    # takes in a v.f. function v(t)(x) which depends on time, and a mesh m
    # also takes in times to evaluate v at
    # returns the persistent barcode for that v(x, t)
    
    P = []
    E = []
    S = []
    
    filename = "pfclmo_v4"
    
    for i, t in enumerate(times):
        print("Processing time " + str(t))
        print("Getting theta...")
        theta_i = CMVF(m, v(t), np.pi/10.0, 0.01)
        
        print("Getting inv...")
        Si = inv(theta_i)
        S += [Si]
        print("Getting pf(cl(Si))...")
        P += [pf(cl(Si), theta_i)]
        print("Getting pf(mo(Si))")
        E += [pf(mo(Si), theta_i)]
        
        # make sure everything is correct subsets
        assert(Si.issubset(P[-1]))
        assert(E[-1].issubset(P[-1]))
    
    print("Intersecting and interleaving the sequences...")
    # compute the intersection between the P and Es
    zagsP = []
    zagsE = []
    zagtimes = []
    for i in range(len(P) - 1):
        zagsP += [P[i] & P[i+1]]
        zagsE += [E[i] & E[i+1]]
        zagtimes += [(times[i] + times[i+1])/2.0]
    
    # interleave the sequences to make it a zigzag persistence, also update the times 
    zigzagP = []
    zigzagE = []
    zigzagtimes = []
    for i in range(len(P) - 1):
        zigzagP += [P[i]]
        zigzagE += [E[i]]
        zigzagtimes += [times[i]]
        zigzagP += [zagsP[i]]
        zigzagE += [zagsE[i]]
        zigzagtimes += [zagtimes[i]]
    zigzagP += [P[-1]]
    zigzagE += [E[-1]]
    zigzagtimes += [times[-1]]
    times = zigzagtimes
    
    # plotting code.  comment out if unneeded
    print("Plotting and saving to file:")
    for i in range(len(zigzagtimes)):
        fig, ax = plt.subplots()
        m.draw2D(ax, rand=False) # draw the mesh
        
        # draw the quiverplot
        xv, yv = np.meshgrid(np.linspace(-7, 7, 15), np.linspace(-7, 7, 15))
        coords = np.stack([xv, yv, np.zeros_like(xv)], axis=2)
        vecs = np.apply_along_axis(v(zigzagtimes[i]), 2, coords)
        plt.quiver(xv, yv, vecs[:, :, 0], vecs[:, :, 1])
        
        for a in zigzagP[i]:
            ax.add_patch(plt.Circle(a.center()[:2], 0.14, color=(0, 0, 1)))
        for b in zigzagE[i]:
            ax.add_patch(plt.Circle(b.center()[:2], 0.12, color=(1, 0.5, 0)))
        if i % 2 == 0:
            for c in S[i//2]:
                ax.add_patch(plt.Circle(c.center()[:2], 0.1, color=(1, 0, 0)))
        
        plt.savefig("./frames/" + filename + "_%02d.png" % i, dpi=300)
        plt.close()
    
    os.chdir("frames")
    subprocess.call([
        'ffmpeg', '-y', '-framerate', '8', '-i', filename + '_%02d.png', '-r', '30', '-pix_fmt', 'yuv420p',
        filename + '.mp4'
    ])
    
    for file_name in glob.glob("*.png"):
        os.remove(file_name)
    
    # assert Ps and Es are a zigzag sequence
    print("Checking zigzagP/E...")
    for i in range(len(zigzagP)):
        if i % 2 == 1:
            assert(zigzagP[i].issubset(zigzagP[i-1]))
            assert(zigzagP[i].issubset(zigzagP[i+1]))
            assert(zigzagE[i].issubset(zigzagE[i-1]))
            assert(zigzagE[i].issubset(zigzagE[i+1]))
    
    # sanity check
    assert(len(zigzagP) == len(P) + len(zagsP))
    assert(len(zigzagE) == len(E) + len(zagsE))
    assert(len(zigzagP) == len(times))
    
    print("Constructing the sequence in Dionysus 2 format...")
    # invent a new vertex which is the dummy vertex
    dummy = len(m.V) + 1
    
    # construct the sequence for use in dionysus' zigzag persistence
    rel_complex = []
    for i in range(len(zigzagP)):
        curr_simplices = set()
        # add all the simplices in P to the current simplices
        for a in zigzagP[i]:
            curr_simplices |= {tuple(a.simplex())}
        
        # for every simplex in E, add that simplex + the dummy vertex to the simplices
        for b in zigzagE[i]:
            assert(tuple(b.simplex()) in curr_simplices)
            curr_simplices |= {tuple(b.simplex() + [dummy])}
            
        rel_complex += [curr_simplices]
        
    all_simplices = set().union(*rel_complex)
    
    # assert that this is a zigzag complex
    print("Checking rel_complex...")
    for i in range(len(rel_complex)):
        assert(len(rel_complex[i]) > 0)
        if i % 2 == 1:
            assert(rel_complex[i].issubset(rel_complex[i-1]))
            assert(rel_complex[i].issubset(rel_complex[i+1]))
            
    # assert that there are simplices
    assert(len(all_simplices) > 0)
    
    # make sure the boundary of every simplex is in the filtration
    to_add = set()
    for s in all_simplices:
        for a in powerset(s):
            if not a is None and not tuple(a) in all_simplices and not tuple(list(a)[::-1]) in all_simplices:
                to_add |= {tuple(sorted(tuple(a)))}
    all_simplices |= to_add
    all_simplices = list(all_simplices)
    
    f = open("simplices.txt", 'w')
    f.write(str(all_simplices))
    f.close()
    
    # finally, get a list of times each simplex was added/removed
    all_simplex_times = []
    for s in all_simplices:
        s_times = []
        for i in range(len(times)):
            t = times[i]
            born = (s in rel_complex[i] and i == 0) or (s in rel_complex[i] and not s in rel_complex[i-1])
            died = (not s in rel_complex[i] and s in rel_complex[i-1])
            if born or died:
                s_times += [t]
        all_simplex_times += [s_times]
        
    f = open("simplextimes.txt", 'w')
    f.write(str(all_simplex_times))
    f.close()
    
    def detail(i,t,d,zz,cells):
        print(i,t,d)
        # if t == 0.0:
        #     for z in zz:
        #         print(z, ' -> ', ' + '.join("%d * (%s)" % (x.element, f[cells[x.index]]) for x in z))
    
    assert(len(all_simplices) == len(all_simplex_times))
    
    print("Creating the filtration...")
    # finally, compute the filtration using Dionysus 2
    f = d.Filtration(all_simplices)
    d.is_simplicial(f, report=True)
    print("Found to be simplicial.  Computing persistence...")
    zz, dgms, cells = d.zigzag_homology_persistence(f, all_simplex_times, callback=detail)
    
    print("Done!  Here's the retrieved diagram:")
    # print something out about the retrieved diagram
    for i, dgm in enumerate(dgms):
        print("Dimension:", i)
        for p in dgm:
            print(p)
            
    # return everything
    return zz, dgms, cells
    
def v1(t):
    # a decent test vector field which has an attracting fixed point until t = 0, then has a repelling fixed point
    d = t + 3.5
    if t < 0:
        def v(x):
            # attracting fixed point
            if np.linalg.norm(x) < d:
                return np.cross(x, np.array([0, 0, 1]))
            else:
                return -x
        return v
    else:
        def v(x):
            # repelling fixed point, some attracting on the outside
            if np.linalg.norm(x) < d:
                return 2*x #+ np.cross(x, np.array([0, 0, 1]))
            else:
                return -x
        return v

def v2(t):
    d = t + 3.5
    def v(x):
        if np.linalg.norm(x) < d:
            return 2*x
        else:
            return -x
    return v

def v3(t):
    d = t + 3.5
    def v(x):
        return -x
    return v

def v4(t):
    d = 2
    def v(x):
        if np.linalg.norm(x) < d:
            return np.cross(x, np.array([0, 0, 1]))
        else:
            return x
    return v

if __name__ == '__main__':
    n = 15
    vert_list = []
    for y in range(n):
        for x in range(n):
            vert_list += [(x - (n-1)/2.0, y - (n-1)/2.0, 0)]
    vert_list = np.vstack(vert_list)
    
    face_list = []
    for i in range(n-1): # row
        for j in range(n-1): # col
            index = i*n + j
            face_list += [(index, index+1, index+n+1)]
            face_list += [(index, index+n+1, index+n)]
    face_list = np.vstack(face_list)
    
    m = hm.Mesh(vert_list, face_list)
    
    fig, ax = plt.subplots()
    # m.draw2D(ax, rand=False) # draw the mesh
    
    # get the persistence diagram and graph it
    dimension_col = {0:(0, 0, 1), 1:(1, 0, 0), 2:(0.5, 1, 0.5)}
    # times = [-2.0, -1.66, -1.33, -1.0, -0.66, -0.33, 0.0, 0.33, 0.66, 1.0, 1.33, 1.66, 2.0, 2.33, 2.66, 3.0]
    # times = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    times = [0.0, 1.0]
    zz, dgms, cells = conley_index_persistence(v4, m, times)
    for i, dgm in enumerate(dgms):
        d.plot.plot_bars(dgm, color=dimension_col[i])
        plt.show()
    
    # define a function on the vertex positions
    # v = lambda x: np.cross(x, np.array([0, 0, 1])) - x
    # v = lambda x: x if np.linalg.norm(x) < 5 else -x
    
    # # draw v
    # xv, yv = np.meshgrid(np.linspace(-7, 7, 15), np.linspace(-7, 7, 15))
    # coords = np.stack([xv, yv, np.zeros_like(xv)], axis=2)
    # vecs = np.apply_along_axis(v, 2, coords)
    # plt.quiver(xv, yv, vecs[:, :, 0], vecs[:, :, 1])
    
    # v_bar = D(m, v, np.pi/10.0, 0.01)
    
    # # draw the arrows
    # for vtx in v_bar.keys():
    #     src = vtx.center()
    #     if not v_bar[vtx] is None:
    #         tgt = v_bar[vtx].center()
    #     ax.arrow(src[0], src[1], tgt[0]-src[0], tgt[1]-src[1], color=(0, 0, 1.0), head_length=0.003*1.5, head_width=0.003*3)
    
    # theta = CMVF(m, v, np.pi/10.0, 0.01, outer_ring=False)
    
    # # get preimage of theta
    # preim = {}
    # for x in theta.keys():
    #     if theta[x] in preim.keys():
    #         preim[theta[x]] += [x]
    #     else:
    #         preim[theta[x]] = [x]
    
    # print("theta 0, 0")
    # cent = list(filter(lambda x: np.linalg.norm(x.center()) <= 0, m.V))[0]
    # print(len(preim[cent]))
    # print(theta[cent].pos[:2])
    
    # choose random element from the image of theta and visualise pf(A)
    # mvA = inv(theta)
    # pfclA = pf(cl(mvA), theta)
    # pfmoA = pf(mo(mvA), theta)
    # for a in pfclA:
    #     ax.add_patch(plt.Circle(a.center()[:2], 0.2, color=(1, 0.5, 0)))
    # for a in pfmoA:
    #     ax.add_patch(plt.Circle(a.center()[:2], 0.2, color=(0.5, 1, 0)))
    # for a in mvA:
    #     ax.add_patch(plt.Circle(a.center()[:2], 0.1, color=(0, 0, 1)))
    
    # choose a random element and visualise mo(A)
    # mvA = set(preim[random.choice(list(preim.keys()))])
    # moA = mo(mvA)
    # for a in mvA:
    #     ax.add_patch(plt.Circle(a.center()[:2], 0.1, color=(0, 0, 1)))
    # for a in moA:
    #     ax.add_patch(plt.Circle(a.center()[:2], 0.1, color=(1, 0.5, 0)))
    
    
    # # plot the preimage of theta
    # for y in preim.keys():
    #     col = (random.random(), random.random(), random.random())
    #     # if len(preim[y]) == 1:
    #     #     h, w = (0.2, 0.2)
    #     #     ax.add_patch(plt.Rectangle(y.center()[:2] - np.array([w/2.0, h/2.0]), w, h, color=(1, 0, 0)))
    #     # if is_critical(set(preim[y])):
    #     #     ax.add_patch(plt.Circle(x.center()[:2], 0.15, color=(1, 0, 0)))
    #     for x in preim[y]:
    #         ax.add_patch(plt.Circle(x.center()[:2], 0.1, color=col))
    
    # for x in inv(theta):
    #     ax.add_patch(plt.Circle(x.center()[:2], 0.15, color=(1, 0, 0)))
    
    # maxlen = 0
    # for y in preim.keys():
    #     if len(preim[y]) > maxlen:
    #         maxlen = len(preim[y])
    # # print(maxlen)
    
    plt.show()