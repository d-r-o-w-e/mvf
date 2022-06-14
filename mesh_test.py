import he_mesh as hm
import numpy as np
import matplotlib.pyplot as plt

m = hm.Mesh()
m.V = [hm.Vertex(np.array([0, 0])),
         hm.Vertex(np.array([1, 0])),
         hm.Vertex(np.array([0, 1]))]
m.E = [hm.Edge(), hm.Edge(), hm.Edge()]
m.F = [hm.Face()]
m.HE = [hm.Halfedge(), hm.Halfedge(), hm.Halfedge(), hm.Halfedge(), hm.Halfedge(), hm.Halfedge()]

# set the vert attributes
m.V[0].he = m.HE[0]
m.V[1].he = m.HE[2]
m.V[2].he = m.HE[4]

# set the edge attributes
m.E[0].he = m.HE[0]
m.E[1].he = m.HE[2]
m.E[2].he = m.HE[4]

# set the face attributes
m.F[0].he = m.HE[0]

# set the halfedge attributes
# 0 done
m.HE[0].twin = m.HE[1]
m.HE[0].next = m.HE[2]
m.HE[0].vert = m.V[0]
m.HE[0].edge = m.E[0]
m.HE[0].face = m.F[0]

# 1 done
m.HE[1].twin = m.HE[0]
m.HE[1].next = m.HE[5]
m.HE[1].vert = m.V[1]
m.HE[1].edge = m.E[0]
m.HE[1].face = m.F[0]

# 2 done
m.HE[2].twin = m.HE[3]
m.HE[2].next = m.HE[4]
m.HE[2].vert = m.V[1]
m.HE[2].edge = m.E[1]
m.HE[2].face = m.F[0]

# 3 done
m.HE[3].twin = m.HE[2]
m.HE[3].next = m.HE[1]
m.HE[3].vert = m.V[2]
m.HE[3].edge = m.E[1]
m.HE[3].face = m.F[0]

# 4 done
m.HE[4].twin = m.HE[5]
m.HE[4].next = m.HE[0]
m.HE[4].vert = m.V[2]
m.HE[4].edge = m.E[2]
m.HE[4].face = m.F[0]

# 5 done
m.HE[5].twin = m.HE[4]
m.HE[5].next = m.HE[3]
m.HE[5].vert = m.V[0]
m.HE[5].edge = m.E[2]
m.HE[5].face = m.F[0]


# edit the mesh to have another vertex, 2 edges, and another face
m.V += [hm.Vertex(np.array([1, 1]))] # 3
m.E += [hm.Edge(), hm.Edge()] # 3, 4
m.F += [hm.Face()] # 1
m.HE += [hm.Halfedge(), hm.Halfedge(), hm.Halfedge(), hm.Halfedge()] # 6, 7, 8, 9

# set new vert attributes
m.V[3].he = m.HE[8]

# set new edge attribs
m.E[3].he = m.HE[6]
m.E[4].he = m.HE[8]

# set new face attribs
m.F[1].he = m.HE[6]

# set new halfedge attribs
# 6 done
m.HE[6].twin = m.HE[7]
m.HE[6].next = m.HE[8]
m.HE[6].vert = m.V[1]
m.HE[6].edge = m.E[3]
m.HE[6].face = m.F[1]

# 7 done
m.HE[7].twin = m.HE[6]
m.HE[7].next = m.HE[1]
m.HE[7].vert = m.V[3]
m.HE[7].edge = m.E[3]
m.HE[7].face = m.F[1]

# 8 done
m.HE[8].twin = m.HE[9]
m.HE[8].next = m.HE[3]
m.HE[8].vert = m.V[3]
m.HE[8].edge = m.E[4]
m.HE[8].face = m.F[1]

# 9
m.HE[9].twin = m.HE[8]
m.HE[9].next = m.HE[7]
m.HE[9].vert = m.V[2]
m.HE[9].edge = m.E[4]
m.HE[9].face = m.F[1]

# fix the other vertices
m.HE[3].next = m.HE[6]
m.HE[5].next = m.HE[9]

n = 10
vert_list = np.array([[np.cos(2*np.pi*i/n), np.sin(2*np.pi*i/n)] for i in range(n)])
face_list = np.array([[0, i, i+1] for i in range(1, n-1)])
m2 = hm.Mesh(vert_list, face_list)


if __name__ == "__main__":
    # m.draw2D(rand=True)
    m2.draw2D(rand=True)
    plt.show()