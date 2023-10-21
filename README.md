# RPS
The purpose of this project is implementing a visibility graph, from a start position to the goal position, in a
2D environment containing obstacles created from polygons by the Rotational Plane Sweep (RPS) algorithm. 
While a brute-force implementation is possible, it has an order $O(n^3)$ but the Rotational Plane Sweep (RPS) algorithm proposed here improves
this to $O\left(n^2 \log(n)\right)$.

You can ﬁnd a description of the RPS algorithm at:
  -     Principles of Robot Motion, H. Choset, et al. Chapter 5.1
  -     https://tanergungor.blogspot.com/2015/04/robot-navigation-rotational-sweep.html

This project is in the studying program of Erasmus Mundus IFRoS students 2023, Autonomous Systems course, Lab 02.

Instructions for contributors about ongoing project. The project contains only one class, named "RPS", implementing the Rotational Plane Sweep (RPS) algorithm by python language.
This class is represented in "RPS.py" file. The main function is at the end of the file. Besides, there are some .csv file, which are the map, where each row [n, x, y] represents
Object number and the x coordinate, y coordinate. How to run the code:
    1.   Run the "RPS.py" file

Result:
    1. Figure 1: The map with the starting point, the goal and the obstacles.
    2. Figure 2: The visibility graph after implementing the Rotational Plane Sweep (RPS) algorithm.
