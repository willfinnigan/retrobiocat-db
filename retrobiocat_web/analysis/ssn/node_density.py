import numpy as np
import time

"""
This is an alternative set of functions for scaling clusters of nodes in an ssn.
Rather than calculating the density, the idea is the find the median (or IQR) of distances between nodes,
and scale based on this.
Ultimately I've not pursued this further, although the code seems to work ok.
"""


def calculate_node_density(pos_dict, around_node, box_distance=100):

    c = pos_dict[around_node]
    x_min, x_max = c[0]-box_distance/2, c[0]+box_distance/2
    y_min, y_max = c[1] - box_distance/2, c[1] + box_distance/2

    nodes_in_box = []
    for node, positions in pos_dict.items():
        x, y = positions[0], positions[1]
        if x >= x_min and x <= x_max and y >= y_min and y <= y_max:
            nodes_in_box.append(node)

    box_size = box_distance * box_distance
    if box_size == 0:
        box_size = 1
    density = len(nodes_in_box) / box_size
    return density


def distance_between_two_nodes(pos_one, pos_two):
    """
    Given two node positions, return the distance between them
    pos_one and pos_two should be in the format (x, y)
    https://stackoverflow.com/questions/1401712/how-can-the-euclidean-distance-be-calculated-with-numpy
    """
    pos_one = np.array(pos_one)
    pos_two = np.array(pos_two)
    return np.linalg.norm(pos_one-pos_two)

def get_closest_node(node_id, node_position, pos_dict):
    """
    Return the closest node to a given node
    """

    current_closest = None
    current_closest_distance = 10**10
    for node, pos in pos_dict.items():
        if node_id != node:
            distance = distance_between_two_nodes(pos, node_position)
            if distance < current_closest_distance:
                current_closest = node
                current_closest_distance = distance

    return current_closest, current_closest_distance

def get_closest_nodes(pos_dict):
    """
    From a dictionary of node positions, calculate the closest distance to another node for every node
    Iterate over pos_dict using get_closest_node to find the closest in each case, saving the results to a dict
    """

    closest = {}
    for node_id, pos in pos_dict.items():
        current_closest, current_closest_distance = get_closest_node(node_id, pos, pos_dict)
        closest[node_id] = current_closest_distance

    return closest

def get_median_distance(pos_dict):
    closest_dict = get_closest_nodes(pos_dict)
    distances = list(closest_dict.values())
    #return np.median(distances)
    return np.percentile(distances, 25)

if __name__ == "__main__":
    t0 = time.time()
    pos_dict = {'test_node_1': (0,0),
                'test_node_2': (1,0),
                'test_node_3': (10,11),
                'test_node_4': (12,13),
                'test_node_5': (5,5)}

    closest = get_closest_nodes(pos_dict)
    t1 = time.time()
    print(closest)
    print(get_median_distance(closest))