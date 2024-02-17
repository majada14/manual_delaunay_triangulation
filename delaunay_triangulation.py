# libraries to import
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import shapely
from scipy.spatial import Delaunay
import numpy as np
import random


# MAIN FUNCTION
def delaunay_triangulation_bowyer_watson(points):
    '''
    this function wraps everything together to compute the Delaunay triangulation given a set of points

    INPUT: a list containing the coordinates of the points
    OUTPUT: a list containing the coordinates of the triangles obtained from the Delaunay algorithm
    '''

    # sort points by the x-coordinate: we want to find the point on the left side.
    points.sort()

    # create a super-triangle large enough to contain all the points
    super_triangle = create_super_triangle(points)
    triangulation = [super_triangle]

    points_added = [] # these structures are used to keep track of what the algorithm is doing
    last_check = [] # these structures are used to keep track of what the algorithm is doing
    for p in points:
        # sort triangles based on the x-coordinate of their leftmost vertex since we proceede in a clockwise direction
        triangulation = clean_triangulation(triangulation)
        triangulation.sort(key=lambda t: min(point[0] for point in t))

        points_added.append(p)

        for triangle in triangulation:
            # sort the vertices by x-coordinate
            list(triangle).sort(key=lambda point: point[0])

            if point_inside_circumcircle(p, triangle):
                ''' here we check if the current triangle contains the current point. 
                If so, we need to remove the triangle and construct new ones '''
                bad_triangles = [triangle]
                edges_to_remove = set()
                triangulation.remove(triangle) # here we remove the triangle from the triangulation
                last_check.append(triangle)

                ''' here we identify the edges of the triangle that we need to remove '''
                for i in range(3):
                    edge = [triangle[i], triangle[(i + 1) % 3]]
                    if not any(edge in t for t in bad_triangles if t != triangle):
                        edges_to_remove.add(tuple(sorted(edge)))

                ''' here we construct the new triangles connecting the point p to the vertices of the triangle in which it is contained '''
                new_triangles = create_new_triangles(p, list(edges_to_remove), triangulation)
                valid_triangles = [t for t in new_triangles if is_triangle(t)] # check the validity of new triangles created
                #valid_triangles = [t for t in new_triangles if is_valid_triangle(t, points_added)]  # check the validity of new triangles created
                triangulation.extend(valid_triangles) # insert the new triangles in the current triangulation

        ''' here we check if new triangles contain some of the points already inserted.
        If so, we identify two adjacent triangles and apply the edge flipping.'''
        wrong_triangles = [t for t in triangulation if not is_valid_triangle(t, points_added)]
        adjacent_triangles = find_adjacent_invalid_triangle_pairs(wrong_triangles) # identify adjacent triangles

        for couple in adjacent_triangles:
            flipped_triangles = handle_invalid_triangles(couple) # flip the triangles
            triangulation.extend(flipped_triangles)  # the flipped triangles are added to the current triangulation

    ''' here we apply some cleaning procedure to the triangulation to ensure that the triangles we are storing are Delaunay triangles. '''
    triangulation = [t for t in triangulation if not any(v in super_triangle for v in t)] # remove the triangles containing a vertex from the original super-triangle
    triangulation = clean_triangulation(triangulation)

    cleaned_triangulation = []
    for t in triangulation:
        to_count = []
        for p in points:
            if point_inside_circumcircle(p, t):
                # if the triangle contains the point
                to_count.append(True)
            else:
                to_count.append(False)

        if True not in to_count:
            # so if the triangle doesn't contain any of the points it's a Delaunay triangle and we should keep it
            cleaned_triangulation.append(t)

    print(f"Delaunay triangulation completed! {len(cleaned_triangulation)} triangles identified.")
    return cleaned_triangulation  # here we return the cleaned Delaunay triangulation

# ------------------------------------------------------------------------------------------------------
# FUNCTIONS TAKING POINTS AS INPUT

def create_super_triangle(points):
    '''
    this function construct a big triangle such that it contains all the points

    INPUT: sorted list of points
    OUTPUT: coordinates of the big triangle
    '''
    min_x = min(p[0] for p in points) # identify the min x in the space
    min_y = min(p[1] for p in points) # identify the min y in the space
    max_x = max(p[0] for p in points) # identify the max x in the space
    max_y = max(p[1] for p in points) # identify the max y in the space

    dx = max_x - min_x
    dy = max_y - min_y
    delta = max(dx, dy)

    p1 = (min_x - 2 * delta, min_y - delta)
    p2 = (min_x + delta, max_y + 2 * delta)
    p3 = (max_x + 2 * delta, min_y - delta)
    return [p1, p2, p3]


def distance_squared(p1, p2):
    '''
    this function compute the euclidean distance between point p1 and point p2

    INPUT:
        - p1: coordinates of point p1
        - p2: coordinates of point p2
    OUTPUT: distance between the two points
    '''
    return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2
# ------------------------------------------------------------------------------------------------------
# FUNCTIONS TAKING TRIANGLES AS INPUT

def circumcenter_of_triangle(triangle):
    '''
    this function identify the circumcenter of the circumcircle obtained from the triangle through the mathematical formula.

    INPUT: list of points of the triangle
    OUTPUT: coordinates of the circumcenter, or None if the points are collinear
    '''
    x1, y1 = triangle[0]
    x2, y2 = triangle[1]
    x3, y3 = triangle[2]
    D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

    if D != 0:
        Ux = ((x1 ** 2 + y1 ** 2) * (y2 - y3) + (x2 ** 2 + y2 ** 2) * (y3 - y1) + (x3 ** 2 + y3 ** 2) * (y1 - y2)) / D
        Uy = ((x1 ** 2 + y1 ** 2) * (x3 - x2) + (x2 ** 2 + y2 ** 2) * (x1 - x3) + (x3 ** 2 + y3 ** 2) * (x2 - x1)) / D
        return (Ux, Uy)
    else:
        # D=0 indicates that the points are collinear and therefore we can't create a triangle
        return None


def find_adjacent_invalid_triangle_pairs(invalid_triangles):
    '''this function identifies adjacent triangles given a set of triangles

    INPUT: invalid triangles: list of non Delaunay triangles
    OUTPUT: list of paired adjacent triangles'''

    pairs = []

    for i in range(len(invalid_triangles)):
        for j in range(i + 1, len(invalid_triangles)):
            common_edge = set(invalid_triangles[i]) & set(invalid_triangles[j])
            if len(common_edge) == 2:  # check if the triangles share an edge
                pairs.append((invalid_triangles[i], invalid_triangles[j]))

    return pairs


def clean_triangulation(triangulation):
    '''this function is used to clean the triangulation from repeated triangles

    INPUT: triangulation: list of triangles
    OUTPUT: new_triangulation: list of triangles'''

    triangulation_set= set() # initialize an empty set
    for triangle in triangulation:
        triangle_object = shapely.Polygon(triangle) # to add elements to a set we need to trandform triangles into objects
        triangulation_set.add(triangle_object)

    if len(triangulation) > len(triangulation_set): # means that some of the triangles are repeated
        new_triangulation = [] # reconstruct the list of triangles with unique triangles
        for triangle in triangulation_set:
            new_triangle = tuple(triangle.exterior.coords[0:3])
            new_triangulation.append(new_triangle)
        return new_triangulation # return the list of triangles without repetitions

    if len(triangulation) == len(triangulation_set): # means that there aren't repeated triangles
        return triangulation # return the original triangulation list


def is_triangle(triangle):
    '''
    Check if the given set of points forms a triangle.

    INPUT: triangle: set of points forming a possible triangle
    OUTPUT: True if the triangle is valid, False otherwise
    '''
    # check if there are three edges and if the three points are not overlapping
    return len(set(triangle)) == 3 and not are_collinear(triangle)


def handle_invalid_triangles(invalid_triangles):
    '''
    This function handles two non Delaunay triangles and transform them into Delaunay triangles through edge flipping

    INPUT: invalid_triangles: list of non Delaunay triangles
    OUTPUT: list of Delaunay triangles (obtained after flipping edge)
    '''

    common_points = set(invalid_triangles[0]) & set(invalid_triangles[1]) # identify the common points between the two triangles
    common_points_list= list(common_points) # turns them into a list
    common_edge = (common_points_list[0], common_points_list[1]) # construct the common edge

    flipped_edge = (tuple(list(set(invalid_triangles[0]) - common_points)[0]), tuple(list(set(invalid_triangles[1]) - common_points))[0]) # construct the flipped edge
    new_triangles = [sorted((flipped_edge[0], flipped_edge[1], common_edge[0])), sorted((flipped_edge[0], flipped_edge[1], common_edge[1]))] # construct the new triangles
    return new_triangles


def are_collinear(triangle):
    '''
    fucntion to check if the points in the given triangle are collinear.

    INPUT: triangle: set of points forming a triangle
    OUTPUT: True if the points are collinear, False otherwise
    '''
    points = list(triangle)  # Convert the set to a list
    x1, y1 = points[0]
    x2, y2 = points[1]
    x3, y3 = points[2]

    return (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) == 0 # to identify collinearity
# ------------------------------------------------------------------------------------------------------
# FUNCTION TAKING MIXED INPUTS

def point_inside_circumcircle(p, triangle):
    '''
    This function checks whether the circumcircle of a triangle contains point p.
    To do so, we compare the distance between the point p and the circumcenter with the radius of the circumcircle.

    INPUT:
        - p: coordinates of a point
        - triangle: set of points constructing a triangle
    OUTPUT: True if the point is inside the circumcircle, False otherwise
    '''

    # if the point is one of the vertices of the triangle we don't consider it
    if p == list(triangle[0]) or p == list(triangle[1]) or p == list(triangle[2]):
        return False

    # here we identify the circumcenter of the circumcircle obtained for the triangle
    circumcenter = circumcenter_of_triangle(triangle)

    # here we check for collinearity between the points
    if circumcenter is None:
        return False

    # if the circumcenter is an integer we continue with the computation
    radius_squared = distance_squared(circumcenter, triangle[0]) # compute the radius of the circumcircle
    distance_squared_to_circumcenter = distance_squared(p, circumcenter) # compute the distance between point p and the circumcenter

    '''consideration:
    if the distance between the circumcenter of the triangle and point p is smaller than the radius of the circumcircle 
    of the same triangle it means that the point is inside the circumcircle of the triangle'''
    return distance_squared_to_circumcenter <= radius_squared


def is_shared_edge(edge, triangles):
    '''
    this function check if there are adjacent triangles

    INPUT:
        - edge: it's the edge we want to check if it is shared between various triangles
        - triangles: list of triangles
    OUTPUT: True if the triangles share the edge, False otherwise
    '''
    # Count how many triangles share the given edge
    count = sum(edge in triangle for triangle in triangles)
    return count == 1


def get_edges_to_flip(p, triangle):
    '''
    this function identify the edge to flip given a rectangle

    INPUT:
        - p: coordinates of a point
        - triangle: list of points constructing the triangle
    OUTPUT: the edge we want to flip
    '''
    edges_to_flip = []
    for i in range(3):
        edge = [triangle[i], triangle[(i + 1) % 3]]
        if is_delaunay_edge(p, edge, triangle):
            edges_to_flip.append(edge)

    return edges_to_flip


def is_delaunay_edge(p, edge, triangle):
    '''
    this function check whether the edges of the current triangulation are delaunay edges or not

    INPUT:
        - p: coordinates of a point
        - edge: line between two points
        - triangle: list of the points composing the triangle
    OUTPUT: True or false if the edge is delaunay or not
    '''

    circumcenter = circumcenter_of_triangle(triangle) # identify the circumcenter of the circumcircle obtained from the triangle
    distance_squared_to_circumcenter = distance_squared(p, circumcenter) # compute the euclidean distance between point p and the circumcenter
    radius_squared = distance_squared(circumcenter, edge[0]) # compute the radius of the circumcircle

    return distance_squared_to_circumcenter <= radius_squared


def create_new_triangles(p, edges_to_flip, triangles):
    '''
    this function construct new triangles given a list of non Delaunay triangles

    INPUT:
        - p: coordinates of a point
        - edges_to_flip: list of the edges that is necessary to flip
        - triangle: list of the non Delaunay triangles
    OUTPUT: list of triangles constructed flipping the edges
    '''
    new_triangles = set()

    for edge in edges_to_flip:
        if not is_shared_edge(edge, triangles):
            new_triangle = [edge[0], edge[1], tuple(p)]
            if is_triangle(new_triangle):
                new_triangles.add(tuple(sorted(new_triangle)))
        else:
            new_triangle1 = list(set(edge) | {tuple(p)})
            if is_triangle(new_triangle1):
                new_triangles.add(tuple(sorted(new_triangle1)))

    return list(sorted(new_triangles))


def is_valid_triangle(triangle, original_points):
    '''
    This function check if any of the original list of points is inside the current triangle

    INPUT:
        - triangle: set of points forming a triangle
        - original_points: list of points added to the set
    OUTPUT: True if the triangle doesn't contain any of the points, False otherwise
    '''

    current_vertices = [list(point) for point in list(triangle)] # here we exclude the vertices of the current triangle
    valid_points = [point for point in original_points if point not in current_vertices]

    if len(valid_points)==0:
        # if any of the original points is inside the triangle the triangle is not valid
        return True

    points_in_triangle = []
    for point in valid_points:
        if point_inside_circumcircle(point, triangle):
            # if the triangle contains a point the triangle is not valid
            points_in_triangle.append(False)
        else:
            # if the triangle doesn't contain a point the triangle is valid
            points_in_triangle.append(True)

    if False in points_in_triangle:
        # if there is at least one point in the triangle, then the triangle is not valid
        return False
    else:
        return len(set(triangle)) == 3 and not are_collinear(triangle)
# ------------------------------------------------------------------------------------------------------
# FUNCTION TO PLOT

def plot_triangles(triangles):
    '''
    function to plot the triangles obtained

    INPUT: triangles: list of triangles
    OUTPUT: plot of the triangles
    '''
    fig, ax = plt.subplots()

    for triangle in triangles:
        polygon = Polygon(triangle, edgecolor='lightblue', facecolor='none')
        ax.add_patch(polygon)

        # Plot vertices in green
        for vertex in triangle:
            ax.plot(vertex[0], vertex[1], 'go', markersize=6)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.title("Delaunay triangulation with manual algorithm")
    plt.savefig(fname="manual_delaunay.png", format="png")
    plt.show()
  

# ------------------------------------------------------------------------------------------------------
''' here we test the algorithm implemented on a set of points and compare it with the results of an external library 
(to check if everything went correctly)'''

random.seed(123) # for reproducibility of the results
# it's important to define three things:

min_coord = [0.0, 0.0] # the minimum coordinates of the plane we use
max_coord = [10.0, 10.0] # the maximum coordinates of the plane we use
num_points = 7 # the number of points

# here we randomly generate the coordinates of the points ww want to analyze
random_points = [[random.uniform(min_coord[0], max_coord[0]), random.uniform(min_coord[1], max_coord[1])] for _ in range(num_points)]

# plot the points
plt.scatter(x=[el[0] for el in random_points],
            y=[el[1] for el in random_points])
plt.title(f"Random set of {len(random_points)} points")
plt.savefig(fname="orignal_points_set.png", format="png")
plt.show()


# now we compute the Delaunay triangulation with our implementation
delaunay = delaunay_triangulation_bowyer_watson(random_points)
plot_triangles(delaunay) # and we plot it

points_array = np.array(random_points)
# now we compute the Delaunay triangulation with the scipy library
delaunay_with_lib = Delaunay(points_array)
plt.triplot(points_array[:,0], points_array[:,1], delaunay_with_lib.simplices) # and we plot it
plt.plot(points_array[:,0], points_array[:,1], 'o')
plt.title("Delaunay triangulation with scipy library")
plt.savefig(fname="scipy_delaunay.png", format="png")
plt.show()
