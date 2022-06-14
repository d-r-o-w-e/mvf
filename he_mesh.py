import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from matplotlib import patches as ptc
import random

class Vertex:
    def __init__(self, pos):
        self.pos = pos
        self.he = None
        self.idx = None
        
    def center(self):
        return self.pos
    
    def simplex(self):
        return [self.idx]

class Edge:
    def __init__(self):
        self.he = None
        
    def center(self):
        return (self.he.vert.pos + self.he.twin.vert.pos)/2.0
    
    def simplex(self):
        return [self.he.vert.idx, self.he.twin.vert.idx]
        
class Face:
    def __init__(self):
        self.he = None
        
    def center(self):
        return (self.he.vert.pos + self.he.next.vert.pos + self.he.next.next.vert.pos)/3.0
    
    def simplex(self):
        return [self.he.vert.idx, self.he.next.vert.idx, self.he.next.next.vert.idx]

class Halfedge:
    def __init__(self):
        self.twin = None
        self.next = None
        self.vert = None
        self.edge = None
        self.face = None
        
class Mesh:
    def __init__(self, vert_list= None, face_list=None):
        # vert list is np array of positions, n x 3
        # face list is np array of indices for each face, n x 3 (0 indexed)
        self.V = []
        self.E = []
        self.F = []
        self.HE = []
        
        if not vert_list is None and not face_list is None:
            # only triangle meshes for now.  sorry
            assert(vert_list.shape[1] == 3)
            assert(face_list.shape[1] == 3)
            
            edge_map = {}
            
            # at least create the vertices
            for i in range(vert_list.shape[0]):
                self.V += [Vertex(vert_list[i, :])]
                self.V[i].idx = i
            
            for i in range(face_list.shape[0]):
                new_face = Face()
                self.F += [new_face]
                for j in range(face_list.shape[1]):
                    new_he = Halfedge()
                    self.HE += [new_he]
                    curr_idx = (int(face_list[i, j]), int(face_list[i, (j+1)%3]))
                    new_he.vert = self.V[curr_idx[0]] # assign vert of he
                    self.V[curr_idx[0]].he = new_he # assign he of vert
                    edge_map[curr_idx] = new_he
                    new_face.he = new_he # assign he of the face
                    new_he.face = new_face

                for j in range(face_list.shape[1]):
                    # set nexts of he
                    curr_idx = (int(face_list[i, j]), int(face_list[i, (j+1)%3]))
                    next_idx = (int(face_list[i, (j+1)%3]), int(face_list[i, (j+2)%3]))
                    edge_map[curr_idx].next = edge_map[next_idx]
                    
                    # set twins of he if possible
                    opp_idx = (curr_idx[1], curr_idx[0])
                    if opp_idx in edge_map.keys():
                        edge_map[curr_idx].twin = edge_map[opp_idx]
                        edge_map[opp_idx].twin = edge_map[curr_idx]
                    
            # fix all the halfedges with no brother.  assumes mesh isn't degenerate
            next_map = {}
            boundary_hes = []
            for he in edge_map.keys():
                if (he[1], he[0]) not in edge_map.keys():
                    new_he = Halfedge()
                    new_he.twin = edge_map[he] # assign new twin
                    edge_map[he].twin = new_he # assign old twin
                    new_he.vert = self.V[he[1]] # assign new vert
                    boundary_hes += [(new_he, he[0], he[1])]
                    self.HE += [new_he]
                    next_map[he[1]] = new_he
            for he, target, source in boundary_hes:
                edge_map[(source, target)] = he # also add to the edge map
                he.next = next_map[target]
                
            # finally, build the edges
            for e in edge_map.keys():
                # only add an edge if it has v1 < v2 or it's on the boundary
                if e[0] < e[1]:
                    new_edge = Edge()
                    self.E += [new_edge]
                    new_edge.he = edge_map[e] # assign he of edge
                    new_edge.he.edge = new_edge
                    new_edge.he.twin.edge = new_edge                                    
        
    def draw2D(self, ax, rand=False):
        # assumes that Z=0 for all elements in the mesh.  does not actually do plt.show()
        # fig, ax = plt.subplots()
        
        # draw all the faces (very light)
        for f in self.F:
            verts = [[f.he.vert.pos[0], f.he.vert.pos[1]], 
                     [f.he.next.vert.pos[0], f.he.next.vert.pos[1]], 
                     [f.he.next.next.vert.pos[0], f.he.next.next.vert.pos[1]]]
            r1 = random.random()*0.2
            r2 = random.random()*0.2
            r3 = random.random()*0.2
            ax.add_patch(ptc.Polygon(verts, closed=True, fill=True, 
                                     color=(1.0-r1, 1.0-r2, 1.0-r3) if rand else (1.0, 1.0, 0.9), zorder=0))
        
        # draw all the edges
        lines = [[(e.he.vert.pos[0], e.he.vert.pos[1]), 
                  (e.he.twin.vert.pos[0], e.he.twin.vert.pos[1])] for e in self.E]
        lc = mc.LineCollection(lines, linewidths=0.5, zorder=1)
        ax.add_collection(lc)
        
        # draw all the points
        x, y = zip(*[(v.pos[0], v.pos[1]) for v in self.V])
        ax.scatter(x, y, s=0.5, c='black', zorder=2)
        
        ax.set_aspect('equal', adjustable='box')
        
        return ax