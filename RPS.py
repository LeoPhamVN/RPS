import numpy as np
import matplotlib.pyplot as plt
import csv

# Class Rotational Plane Sweep
class RPS:
    """
    RPS class implements the Rotation Plane Sweep algorithm to compute the visibility 
    of the starting point with respect to all the vertices and goal point.


   - vertices: array of size N x 3, where each row [n, x, y] represents
     Object number and the x coordinate, y coordinate. 
     The first row represents the START vertex and the last one the GOAL.
   - edges: array of size M x 2, where each row [vi ve] represents the
     initial vertex and the final vertex of the edge. M is the number of
     edges, which is not pre-defined, given that it depends on the problem.
    """

    """
    Public Methods
    """
    def __init__(self, vertices:np.array(float)) -> None:
        """
        Init function
        ...
          
        ...

        Parameters
        -------------
        Private variables:

        __vertices  : an array [N x 3]
                    Array contains the number of the polygon the point is belongs to and the position of points in the map.
                    The first row is the starting point and the last row is the goal.
                    N: The number of point.

        __S         : an array
                    Sorted list of edges: S
        
        __v         : an array [1 x 3]
                    Vertex v: start point

        __vi        : an array [1 x 3]
                    Vertex vi: whose visibility will be tested

        __E         : a list
                    Edge list size M x 2, where M is the number of edges. Each row represents one edge, having the start and the end vertices index of the array of vertices.
        
        __E_dst     : an array
                    Distance
        
        __edges     : a list
                    List edges in the visibility list. Each cell contains 2 elements, indexs of the starting point and the end point.
        
        __edges_idx : int
                    Index of the next edge to be inserted.

        __A         : an array [1 x N-1]
                    Array contains the indices of the elements of softed alpha.
        Returns
        -------------

        """
        self.__vertices     = vertices              
        self.__S            = np.array(float)       
        self.__v            = np.array(float)       
        self.__vi           = np.array(float)       
        self.__E            = []                    
        self.__E_dst        = np.array(float)       
        self.__edges_idx    = 1                     
        self.__edges        = []                   
        self.__A            = np.array(int)         

    def get_edge(self):
        """
        A public function used to get visibility edge list
        ...
            
        ...

        Parameters
        -------------
        
        Returns
        -------------
            __edges: list
                The edge list
        """
        return self.__edges

    def compute(self):
        """
        The main function
        ...
            
        ...

        Parameters
        -------------
        self : 
            
        Returns
        -------------
            True
        """
        # Number of vertices
        N = self.__vertices.shape[0]

        # List of edges
        self.__calculate_edges()
        
        # Plot the initial environment
        plt.figure()
        self.__plot_environment()
        
        # Test whether the environment has no obstacles
        if len(self.__E) == 0:
            self.__edges = np.array([[1, 2]])
            self.plot_result()
            return True
        
        # Iterate through all the vertices to determine the visible vertices
        # from each vertex
        for i in range(N):
            # Vertex v: start point
            self.__v = self.__vertices[i]

            # Subset of vertex except the start point
            subset = np.concatenate((self.__vertices[:i], self.__vertices[i+1:]), axis=0)

            # Angle from the horizontal axis to the line segment vv_i sorted in
            # incresing order
            self.__calculate_alpha(subset)

            # Sorted list of edges that intersects the horizontal line
            # emanating from v
            self.__intersects_line()
            
            # Evaluate each vertex
            for j in range(N - 1):

                # Determine the index of the vertex in the initial array
                vertex_nr = self.__A[j] if self.__A[j] < i else self.__A[j] + 1

                # Vertex whose visibility will be tested
                self.__vi = subset[self.__A[j]]

                # Add the edge to the visible list, in case is visible
                self.__add_edge(i, vertex_nr)
                
                # Determine the edges indexes where vi is the start edge and the end edge
                start_edge = np.where(self.__E[:, 0] == vertex_nr)[0]
                end_edge = np.where(self.__E[:, 1] == vertex_nr)[0]
                if len(start_edge) and len(end_edge):
                    start_edge = int(start_edge.item())
                    end_edge = int(end_edge.item())
                    # Find the edges that should be either deleted or inserted
                    insert_edges, delete_edges = self.__find_edges(start_edge, end_edge)
                    
                    # If vi is in the begining of an edge that is not in S
                    if len(insert_edges):
                        self.__insert_edge(insert_edges)
                    
                    # If vi is in the end of an edge in S
                    if len(delete_edges):
                        self.__delete_edge(delete_edges)
        
        # Delete internal edges and add polygon lines
        self.__clear_edges()

        # Plots the resulting edges
        self.__plot_result()

        # Paint in blue the polygons
        self.__plot_environment()

    """
    Private Methods
    """
    def __add_edge(self, i:int, vertex_nr:np.array(int)):
        """
        ADD_EDGE Add visible edges to the vector edges and update the value of the index
        ...
        
        ...
        Parameters
        ---------------
        vertex_nr   :


        Returns
        ---------------

        """

        if len(self.__S):
            # Test whether the vertex is visible
            if self.__is_visible():
                # Add indexes [v, vi] to the visibility graph
                self.__edges.append([i, vertex_nr]) 
                self.__edges_idx += 1
        else:
            # If S is empty, add index to the visibility graph
            self.__edges.append([i, vertex_nr]) 
            self.__edges_idx += 1

    def __find_edges(self, start_idx:int, end_idx:int):
        """
        FIND_EDGES Determines the initial and the end edges
        """
        insert_edges = []
        delete_edges = []

        proj_start = 0
        proj_end = 0

        line_v_vi_hmg = self.__homogeneous_coordinates([self.__v[1], self.__v[2], self.__vi[1], self.__vi[2]])
        line_v_vi_hmg /= np.linalg.norm(line_v_vi_hmg[0:2])

        vertex_start = np.concatenate((self.__vertices[self.__E[start_idx, 1], 1:3], [1]))
        vertex_end = np.concatenate((self.__vertices[self.__E[end_idx, 0], 1:3], [1]))

        if start_idx is not None:
            proj_start = np.dot(line_v_vi_hmg, vertex_start)
        
        if end_idx is not None:
            proj_end = np.dot(line_v_vi_hmg, vertex_end)

        insert_idx = 0
        delete_idx = 0

        if proj_start > 0:
            insert_edges.append(start_idx)
            insert_idx += 1
        if proj_start < 0:
            delete_edges.append(start_idx)
            delete_idx += 1
        
        if proj_end > 0:
            insert_edges.append(end_idx)
            insert_idx += 1
        if proj_end < 0:
            delete_edges.append(end_idx)
            delete_idx += 1

        return np.array(insert_edges), np.array(delete_edges)

    def __calculate_edges(self):
        """
        Transforms the initial set of vertices into a data
        structure that represents the edges of the polygons.
        
        __E: Array of size M x 2, where M is the number of edges. Each row
        represents one edge, having the start and the end vertices index of the
        array of vertices.
        """

        N = self.__vertices.shape[0]
        # Initialization of the variables
        self.__E = []
        edge_idx = 0

        for i in range(1, N - 1):

            # Number of the current polygon
            object_nr = self.__vertices[i, 0]

            # Check previous vertex to find the initial vertex of the object
            if self.__vertices[i - 1, 0] != object_nr:
                # Initial vertex index of this object
                init_vertex_idx = i

            # Check next vertex to find an edge
            if self.__vertices[i + 1, 0] == object_nr:
                # the next vertex belongs to the same object, so there is an edge
                self.__E.append([i, i + 1])
            else:
                # the current vertex is the last one of the object, so there is
                # an edge between this vertex and the initial vertex of the
                # object
                self.__E.append([i, init_vertex_idx])

            edge_idx += 1

        # Size M x 2
        self.__E = np.array(self.__E)

    def __calculate_alpha(self, subset:np.array(float)):
        """
        Calculates the angle from the horizontal axis, to the line
        segment defined by each vertex of the vertices array and v. 
        
        Parameters
        -------------
        - subset: array of size N x 3, where each row [n, x, y] represents
             object number and the x coordinate, y coordinate. 
        """

        # Number of vertices 
        N = subset.shape[0]

        alpha = np.zeros(N)

        # Calculates the angle from the horizontal axis
        for i in range(N):
            x = subset[i, 1] - self.__v[1]
            y = subset[i, 2] - self.__v[2]
            # Compute the angle in the interval [-pi, pi]
            alpha[i] = np.arctan2(y, x)

        # Convert angles to the interval [0, 2*pi]
        alpha = np.mod(alpha, 2 * np.pi)

        # Sort the elements of alpha and return the indices
        # size 1 x N
        self.__A = np.argsort(alpha)

    def __intersects_line(self):
        """
        Select the edges that intersects the horizontal half-line emanating from v
        Initialise the S list with the new vertex.
        """
        # Initialisation of the vectors
        E_dst = []
        S = []

        # Number of edges
        N = self.__E.shape[0]

        # Maximum x-coordinate
        max_idx     = np.argmax(self.__vertices[:, 1])

        # Horizontal half-line emanating from v
        line_v      = [self.__v[1], self.__v[2], self.__vertices[max_idx, 1], self.__v[2]]
        
        # Index of the S array
        s_idx = 0
        
        # Determines the edges that intersects the line
        for i in range(N):

            # Define edge as a line in the form [x1 y1 x2 y2]
            line = np.concatenate((self.__vertices[self.__E[i, 0], 1:3], self.__vertices[self.__E[i, 1], 1:3]))
            
            # determine wheter lines intersects or not and compute the distance
            # to the initial vertex
            intersect, dst = self.__is_intersected(line_v, line)
            
            # Define whether the lines intersects or not
            if intersect:
                # determine if the lines are different segments of the same
                #  line. If this is the case, the edge is not occluding any
                #  line, so it is not included on the S list
                if np.linalg.norm((self.__homogeneous_coordinates(line_v) - self.__homogeneous_coordinates(line))) > 0.00001:
                    # Add edge index to the list
                    S.append(i)
                    # Assign the distance
                    E_dst.append(dst)
                    s_idx += 1

        # Sort the elements of E_dst and return the indexes
        sorted_idx = np.argsort(E_dst)
        self.__E_dst = np.array(E_dst)[sorted_idx]

        # Sort S according to the sorted indexes
        S = np.array(S)[sorted_idx]
        self.__S = S.astype(int)

    def __is_intersected(self, line_1, line_2):
        """
        IS_INTERSECTED Determines whether two lines intersects or not
        - line_1: vector of size 1 x 4, in the way [x1 y1 x2 y2], where x1 and
            y1 is the initial point of the line and x2 y2 the end point
        - line_2: vector of size 1 x 4, in the way [x1 y1 x2 y2], where x1 and
        y1 is the initial point of the line and x2 y2 the end point
        """

        # Initialise the value of the distance
        dst = 0

        # Transform line 1 to homogeneous coordinates
        line_hmg_1 = self.__homogeneous_coordinates(line_1)

        # Transform line 2 to slope-intersect form
        line_hmg_2 = self.__homogeneous_coordinates(line_2)

        # Intersection point
        intersect_pt = np.cross(line_hmg_1, line_hmg_2)

        if intersect_pt[2] == 0:
            intersect = False
        else:
            # Normalize the vector, so the third component is one
            intersect_pt /= intersect_pt[2]

            # X-coordinate of the intersection
            x = intersect_pt[0]

            # sort the x values of each line
            x_line_1 = np.sort([line_1[0], line_1[2]])
            x_line_2 = np.sort([line_2[0], line_2[2]])

            # check if the intersection is on an edge. In this case, there is no occlusion
            if (self.__is_an_edge(intersect_pt, line_1) or
                self.__is_an_edge(intersect_pt, line_2)):
                intersect = False
            else:
                # X-coordinate is on the lines
                if (x >= x_line_1[0] and x <= x_line_1[1] and
                    x >= x_line_2[0] and x <= x_line_2[1]):
                    intersect = True
                    # Euclidean distance
                    dst = np.linalg.norm(intersect_pt[0:2] - line_1[0:2])
                else:
                    intersect = False

        return intersect, dst        

    def __is_an_edge(_, intersect_pt, line):
        """
        IS_AN_EDGE Determines whether the point defined by intersect_pt lies in
        one of the vertex of and edge.
        
        The difference between both coordinates is computed and evaluated agains
        a threshold.
        """
        diff_x_line_1 = abs(intersect_pt[0] - line[0])
        diff_y_line_1 = abs(intersect_pt[1] - line[1])
        diff_x_line_2 = abs(intersect_pt[0] - line[2])
        diff_y_line_2 = abs(intersect_pt[1] - line[3])

        if ((diff_x_line_1 < 0.0001 and diff_y_line_1 < 0.0001) or
            (diff_x_line_2 < 0.0001 and diff_y_line_2 < 0.0001)):
            edge = True
        else:
            edge = False

        return edge

    def __homogeneous_coordinates(_, line):
        """
        HOMOGENEOUS_COORDINATES Transforms a line having the starting and ending
        points to homogeneous coordinates
        
        - line: vector of size 1 x 4, in the way [x1 y1 x2 y2], where x1 and
            y1 is the initial point of the line and x2 y2 the end point
        """
        a = line[1] - line[3]
        b = line[2] - line[0]
        c = (line[0] * line[3]) - (line[2] * line[1])

        hmg = np.array([a, b, c])
        return hmg

    def __is_visible(self):
        """
        S_VISIBLE Determines whether a vertex is visible or not
        """
        # Distance from v to vi
        dst = np.linalg.norm(self.__v[1:3] - self.__vi[1:3], ord=2)

        # Number of edges in S
        N = len(self.__S)

        # Initialize the visible variable
        visible = True

        # Define the line that goes from v to the evaluated vertex, vi
        line_v_vi = [self.__v[1], self.__v[2], self.__vi[1], self.__vi[2]]

        # Determine if the line v_vi intersects with any edge that is closer to
        # the candidate vertex
        S_idx = 0
        while S_idx < N:
            E_idx = self.__S[S_idx]

            line_e = np.concatenate((self.__vertices[self.__E[E_idx, 0], 1:3], self.__vertices[self.__E[E_idx, 1], 1:3]))

            intersect, _ = self.__is_intersected(line_v_vi, line_e)

            if intersect:
                visible = False
                break

            S_idx += 1

        return visible
    
    def __insert_edge(self, insert_edges):
        """
        INSERT_EDGES Summary of this function goes here
        Detailed explanation goes here
        """

        # Number of edges to be inserted
        N = insert_edges.shape[0]
        for i in range(N):
            # The edge is not in S
            if insert_edges[i] not in self.__S:

                # Calculate the distance from the origin to the vertex of the edge
                dst = np.linalg.norm(self.__v[1:3] - self.__vi[1:3], 2)
                
                # Find the indices of elements smaller and bigger than dst in E_dst
                smaller_idx = np.where(self.__E_dst <= dst)
                bigger_idx = np.where(self.__E_dst > dst)
                
                # Insert the edge into S at the appropriate position
                self.__S = np.concatenate((self.__S[smaller_idx], [int(insert_edges[i])], self.__S[bigger_idx]))
                # Insert the distance into E_dst at the appropriate position
                self.__E_dst = np.concatenate((self.__E_dst[smaller_idx], [dst], self.__E_dst[bigger_idx]))

    def __delete_edge(self, edges_dlt):
        """
        INSERT_EDGES Summary of this function goes here
        Detailed explanation goes here
        """
        # Number of edges to be deleted
        N = edges_dlt.shape[0]
        S_size = self.__S.shape[0]
        
        for i in range(N):
            # Check if the edge is in S
            if edges_dlt[i] in self.__S:
                idx_S = np.where(self.__S == edges_dlt[i])
                self.__S = np.delete(self.__S, idx_S)
                self.__E_dst = np.delete(self.__E_dst, idx_S)
        
    def __clear_edges(self):
        """
        UNTITLED6 Summary of this function goes here
        Detailed explanation goes here
        """
        edge_idx = 0
        i = 0

        # Delete edges that belongs to the same polygon
        while i < len(self.__edges):
            if self.__vertices[self.__edges[i][0], 0] == self.__vertices[self.__edges[i][1], 0]:
                self.__edges = np.delete(self.__edges, i, axis=0)
            else:
                edge_idx += 1
                i += 1

        self.__edges = np.concatenate((self.__edges, self.__E), axis=0)
        self.__edges = self.__edges[np.lexsort(np.fliplr(self.__edges).T)]

    def __plot_environment(self):
        """
        Summary of this function goes here detailed explanation goes here
        """
        # Start
        plt.plot(self.__vertices[0, 1], self.__vertices[0, 2], 'ro', markersize=8, markerfacecolor='red')

        # Goal
        plt.plot(self.__vertices[-1, 1], self.__vertices[-1, 2], 'go', markersize=8, markerfacecolor='green')

        # Plot the edges
        for i in range(len(self.__E)):
            plt.plot([self.__vertices[self.__E[i, 0], 1], self.__vertices[self.__E[i, 1], 1]],
                    [self.__vertices[self.__E[i, 0], 2], self.__vertices[self.__E[i, 1], 2]], 'bo', markersize=8, markerfacecolor='blue')
            plt.plot([self.__vertices[self.__E[i, 0], 1], self.__vertices[self.__E[i, 1], 1]],
                    [self.__vertices[self.__E[i, 0], 2], self.__vertices[self.__E[i, 1], 2]], 'b--')

        plt.show()

    def __plot_result(self):
        """
        Summary of this function goes here detailed explanation goes here
        """
        for n in range(self.__edges.shape[0]):
            plt.plot([self.__vertices[self.__edges[n, 0], 1], self.__vertices[self.__edges[n, 1], 1]],
                    [self.__vertices[self.__edges[n, 0], 2], self.__vertices[self.__edges[n, 1], 2]], '-r')


# Main function    
if __name__ == "__main__" :
    # opening the CSV file
    with open('env_4.csv', mode ='r') as file:       
        # reading the CSV file
        csvFile = csv.reader(file, delimiter=',')
        header = []
        header = next(csvFile)
        csvFile = list(csvFile)
        # Get the vertices from the CSV file
        vertices = np.array(csvFile).astype(float)

    # Init object of RPS class
    RPS_object = RPS(vertices)
    # Run main function of RPS class
    RPS_object.compute()
    # Get the list of visibility
    edge = RPS_object.get_edge()