#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 08:24:29 2022

@author: mikkonea
"""

#TODO Tee tämä tuonne Rayscaan. Myös testaa että miten Rayscassa
# voisi interp1d:n sijasta käyttää vastaavaa NN-rakennetta.

import numpy as np

N = 30

data = np.arange(N*N*N).reshape((N,N,N))

dims = [0.2,0.2,0.2]

coords = np.zeros((N*N*N,3))
origin = np.array((0.0, 0.0, 0.0))

linidx = 0
for i in range(N):
    for j in range(N):
        for k in range(N):
            
            coords[linidx,:] = np.array([i * dims[0] + origin[0],
                                         j * dims[1] + origin[1],
                                         k * dims[2] + origin[2]])
            linidx = linidx + 1

def direct_nn(point,coords):
    return np.argmin(np.linalg.norm(coords - point,axis=1))

class NearNeigh_Tree:
    def __init__(self):        
        # The point mask contains all the points in this tree level before
        # making the split
        self.point_mask = None
        
        # The plane splits the remaining medium points in two groups
        # Plane has a point and a normal
        self.plane = None
        
        # Positive/negative child is the nn_tree containing the points
        # on the positive/negative side of the plane (according to the 
        # plane normal).
        self.positive_child = None
        self.negative_child = None

    def get_point_sign(self,point):
        # Return True if point is in the positive side.
        i = self.plane['normal_ind']
        return point[i] >= self.plane['point'][i]
    
    def get_neighbourhood_mask(self,point):
        if self.positive_child == None or self.negative_child == None:
            return self.point_mask
        elif self.get_point_sign(point):
            return self.positive_child.get_neighbourhood_mask(point)
        else:
            return self.negative_child.get_neighbourhood_mask(point)
            
            
def generate_plane(coords):
    mid_point = np.mean(coords,axis=0)
    best_normal = 0
    split_score = np.inf
    for i in range(mid_point.size):
        normal = np.zeros_like(mid_point)
        normal[i] = 1.0
        split_points = coords[:,i] > mid_point[i]
        current_score = np.abs(coords.shape[0] // 2 - np.where(split_points)[0].size)
        if current_score < split_score:
            best_normal = normal
            split_score = current_score
    return {'point': mid_point, 'normal' : best_normal, 
            'normal_ind' : np.where(best_normal)[0][0]}

def generate_masks(coords,initial_mask,plane):
    i = plane['normal_ind']
    m_plus = coords[:,i] >= plane['point'][i]
    m_minus = coords[:,i] < plane['point'][i]
    return (np.logical_and(m_plus,initial_mask),
            np.logical_and(m_minus,initial_mask))
"""
mp = np.array(np.ones((coords.shape[0])),dtype=bool)
for i in range(8):
    p = generate_plane(coords[mp])
    mp,mm = generate_masks(coords,mp,p)
    print(np.where(mp)[0].size)
    print(p)
"""
def populate(coords,tree):
    low_limit = 16 # arbitrary at this point
    if np.where(tree.point_mask)[0].size < low_limit:
        return
    plane = generate_plane(coords[tree.point_mask])
    tree.plane = plane
    mask_p, mask_n = generate_masks(coords, tree.point_mask, plane)
    tree.positive_child = NearNeigh_Tree()
    tree.positive_child.point_mask = mask_p
    tree.negative_child = NearNeigh_Tree()
    tree.negative_child.point_mask = mask_n
    populate(coords,tree.positive_child)
    populate(coords,tree.negative_child)

def create_nn_tree(coords):
    new_tree = NearNeigh_Tree()
    new_tree.point_mask = np.array(np.ones((coords.shape[0])),dtype=bool)
    new_tree.plane = generate_plane(coords)
    populate(coords,new_tree)
    return new_tree

def tree_nn(point,coords,tree):
    mask = tree.get_neighbourhood_mask(point)
    idx_in_mask = np.argmin(np.linalg.norm((coords[mask,:] - point),axis=1))
    return np.where(mask)[0][idx_in_mask]

point = np.array([0.2, 0.0, 0.0])
import time

def get_random_point():
    return np.array([np.random.random() * N * dims[0] + origin[0],
                     np.random.random() * N * dims[1] + origin[1],
                     np.random.random() * N * dims[2] + origin[2]
                     ])
timing_mode = True
n_points = 10000
if timing_mode:

    start_build_tree = time.time()
    tree = create_nn_tree(coords)
    end_build_tree = time.time()
    for i in range(n_points):
        rand_point = get_random_point()
        idx_tree = tree_nn(rand_point,coords,tree)
    end_search_tree = time.time()
    for i in range(n_points):
        rand_point = get_random_point()
        idx_norm = direct_nn(rand_point,coords)
    end_search_norm = time.time()

print("Building tree time: %f" % (end_build_tree - start_build_tree))
print("Search time tree: %f" % (end_search_tree - end_build_tree))
print("Search time normal: %f" % (end_search_norm - end_search_tree))

validity_mode = True
err_list_tree = []
err_list_norm = []
if validity_mode:
    for i in range(n_points):
        rand_point = get_random_point()
        idx_tree = tree_nn(rand_point,coords,tree)
        idx_norm = direct_nn(rand_point,coords)
        err_list_tree.append(np.linalg.norm(coords[idx_tree,:] - rand_point))
        err_list_norm.append(np.linalg.norm(coords[idx_norm,:] - rand_point))
        #if (idx_norm != idx_tree):
        #    print(rand_point,idx_norm,idx_tree)
        #    #TODO: laske virhe!!!
import matplotlib.pyplot as plt
#print(np.mean(err_list))

plt.figure()
plt.hist(err_list_tree,100,alpha=0.5)
plt.hist(err_list_norm,100,alpha=0.5)
plt.show()

if False:
    import pickle5 as pickle
    with open('../../plume_direct.pkl','rb') as f:
        direct_plume = pickle.load(f)
    with open('../../plume_rotated.pkl','rb') as f:
        rotated_plume = pickle.load(f)
    
    import matplotlib.pyplot as plt
    plt.imshow(direct_plume)
    plt.show()
    plt.imshow(rotated_plume)
    plt.show()
    