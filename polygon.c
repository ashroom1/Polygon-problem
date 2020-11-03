#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define INF INT32_MAX
//#define Test 1

/*
 * Check strict math settings
 * We rely on the assumption that we store the vertices in an order where there is an edge from vertex i to i + 1
 * If slope is infinity c is in the form x - c = 0
 * Make sure to cycle from the last vertex to vertex 0
 */

_Bool point_in_polygon(double x, double y);

struct Line {
    double m;
    double c;
    int x1;
    int y1;
    int x2;
    int y2;
}*edge_equations;

struct Line_Double {
    double m;
    double c;
    double x1;
    double y1;
    double x2;
    double y2;
};

struct intersept{
    double x;
    double y;
} *intersepts;

int **vertices;
int num_vertices;
int num_edges;


double r(double x) {
    if (!((int)fabs(x))) {
        if ((x > 0 && x - 0.0000001 < 0) || (x < 0 && x + 0.0000001 > 0))
            return round(x);
    }
    if(((int) (fabs(x) + 0.0000001) > (int) fabs(x)) || ((int) (fabs(x) - 0.0000001) < (int) fabs(x))) {
        return round(x);
    }
    return x;
}

struct Line FindMC(int x1, int y1, int x2, int y2) {

    struct Line line = { .m = (x1 == x2 ? INF : r(((1.0 * (y2 - y1)) / (x2 - x1)))), .c = r(((x1 == x2) ? (0 - x1) : (y1 - (x1 * (1.0 * (y2 - y1) / (x2 - x1))))))};
    return line;
}

void make_edge_equation(int pos_in_array, int x1, int y1, int x2, int y2) {

    if (x1 == x2) {
        edge_equations[pos_in_array].m = INF;
        edge_equations[pos_in_array].c = 0 - x1;
        edge_equations[pos_in_array].x1 = x1;
        edge_equations[pos_in_array].y1 = y1;
        edge_equations[pos_in_array].x2 = x2;
        edge_equations[pos_in_array].y2 = y2;
        return;
    }
    struct Line line = FindMC(x1, y1, x2, y2);
    edge_equations[pos_in_array].m = line.m;
    edge_equations[pos_in_array].c = line.c;
    edge_equations[pos_in_array].x1 = x1;
    edge_equations[pos_in_array].y1 = y1;
    edge_equations[pos_in_array].x2 = x2;
    edge_equations[pos_in_array].y2 = y2;
}

double lines_intersection[2];
double intersection_vertex_number;

_Bool on_same_side(int const point1[], int const point2[], struct Line current_line) {
    double value1 = current_line.m == INF ? r(point1[0] + current_line.c) : r(point1[1] - current_line.m * point1[0] - current_line.c);
    double value2 = current_line.m == INF ? r(point2[0] + current_line.c) : r(point2[1] - current_line.m * point2[0] - current_line.c);
    return ((value1 > 0 && value2 > 0) || (value1 < 0 && value2 < 0)) ? 1 : 0;
}

int is_vertex(double x, double y) {

    for (int i = 0; i < num_vertices; ++i)
        if((double)vertices[i][0] == x && (double)vertices[i][1] == y)
            return i;

    return -1;
}

_Bool consider_point(struct intersept in, int x1, int y1, int x2, int y2, struct Line current_line) {

    double slope = current_line.m;
    int index;
    if ((index = is_vertex(in.x, in.y)) > -1) {
        if (in.x == x1 && in.y == y1) {
            if (!point_in_polygon((1.0 * ((x2 * 1) - (x1 * 1000000))) / (1 - 1000000),
                                  (1.0 * ((y2 * 1) - (y1 * 1000000))) / (1 - 1000000))) {
                return 1;
            } else
                return 0;
        } else if (in.x == x2 && in.y == y2) {
            if (!point_in_polygon((1.0 * ((x1 * 1) - (x2 * 1000000))) / (1 - 1000000),
                                  (1.0 * ((y1 * 1) - (y2 * 1000000))) / (1 - 1000000))) {
                return 1;
            } else
                return 0;
        }

        if (edge_equations[index].m == slope /*|| edge_equations[(index -1) % num_edges].m == slope */||
            edge_equations[(index - 1) % num_edges].m == slope) {
            int index_overlap = -1000;      //Initialized to this to try and make it an invalid array index
            if (edge_equations[index].m == slope) {
                index_overlap = index;
                if (sqrt(pow(x2 - in.x, 2) + pow(y2 - in.y, 2)) < sqrt(pow(x2 - (edge_equations[index].x1 == in.x &&
                                                                                 edge_equations[index].y1 == in.y
                                                                                 ? edge_equations[index].x2
                                                                                 : edge_equations[index].x1), 2) +
                                                                       pow(y2 - (edge_equations[index].x1 == in.x &&
                                                                                 edge_equations[index].y1 == in.y
                                                                                 ? edge_equations[index].y2
                                                                                 : edge_equations[index].y1), 2)))
                    return 0;
            } /*else if (edge_equations[(index + 1) % num_edges].m == slope) {
                index_overlap = (index + 1) % num_edges;
                if (sqrt(pow(x2 - in.x, 2) + pow(y2 - in.y, 2)) < sqrt(pow(x2 - (edge_equations[(index + 1) %
                                                                                                num_edges].x1 == in.x &&
                                                                                 edge_equations[(index + 1) %
                                                                                                num_edges].y1 == in.y
                                                                                 ? edge_equations[(index + 1) %
                                                                                                  num_edges].x2
                                                                                 : edge_equations[(index + 1) %
                                                                                                  num_edges].x1), 2) +
                                                                       pow(y2 - (edge_equations[(index + 1) %
                                                                                                num_edges].x1 == in.x &&
                                                                                 edge_equations[(index + 1) %
                                                                                                num_edges].y1 == in.y
                                                                                 ? edge_equations[(index + 1) %
                                                                                                  num_edges].y2
                                                                                 : edge_equations[(index + 1) %
                                                                                                  num_edges].y1), 2)))
                return 0;
         }*/else if (edge_equations[(index - 1) % num_edges].m == slope) {
                index_overlap = (index - 1) % num_edges;
                if (sqrt(pow(x2 - in.x, 2) + pow(y2 - in.y, 2)) < sqrt(pow(x2 - (edge_equations[(index - 1) %
                                                                                                num_edges].x1 == in.x &&
                                                                                 edge_equations[(index - 1) %
                                                                                                num_edges].y1 == in.y
                                                                                 ? edge_equations[(index - 1) %
                                                                                                  num_edges].x2
                                                                                 : edge_equations[(index - 1) %
                                                                                                  num_edges].x1), 2) +
                                                                       pow(y2 - (edge_equations[(index - 1) %
                                                                                                num_edges].x1 == in.x &&
                                                                                 edge_equations[(index - 1) %
                                                                                                num_edges].y1 == in.y
                                                                                 ? edge_equations[(index - 1) %
                                                                                                  num_edges].y2
                                                                                 : edge_equations[(index - 1) %
                                                                                                  num_edges].y1), 2)))
                    return 0;
            }
            int point1[2], point2[2];
            point1[0] = edge_equations[(index_overlap - 1) % num_edges].x1;
            point1[1] = edge_equations[(index_overlap - 1) % num_edges].y1;
            point2[0] = edge_equations[(index_overlap + 1) % num_edges].x2;
            point2[1] = edge_equations[(index_overlap + 1) % num_edges].y2;

            if (on_same_side(point1, point2, current_line))
                return 0;
            else
                return 1;
        } else {
            int point1[] = {edge_equations[(index - 1) % num_edges].x1, edge_equations[(index - 1) % num_edges].y1};
            int point2[] = {edge_equations[index].x2, edge_equations[index].y2};
            if (on_same_side(point1, point2, current_line))
                return 0;
            else
                return 1;
        }
    }

    return 1;
}


double find_length_using_intersepts(int intersept_counter, int x1, int y1, int x2, int y2, struct Line current_line) {

    double sideA = DBL_MAX, sideB = DBL_MAX;

#ifdef Test

    for (int j = 0; j < intersept_counter; ++j)
        printf("%lf %lf\n", intersepts[j].x, intersepts[j].y);

#endif

    for (int i = 0; i < intersept_counter; ++i) {

        if (consider_point(intersepts[i], x1, y1, x2, y2, current_line)) {
                double d1 = sqrt(pow((intersepts[i].y - y1), 2) + pow((intersepts[i].x - x1), 2));
                double d2 = sqrt(pow((intersepts[i].y - y2), 2) + pow((intersepts[i].x - x2), 2));

                if (d1 == 0.0)
                    sideA = 0;
                if (d2 == 0.0)
                    sideB = 0;

                //if (d1 > d2) {
                if (d1 < d2) {
                    if (sideA > d1)
                        sideA = d1;
                }
                //else if (d2 > d1) {
                else if (d2 < d1) {
                    if (sideB > d2)
                        sideB = d2;
                }
        }
    }

    return sideA + sideB + sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}


int find_line_intersection(struct Line_Double l1, struct Line l2) //0 if intersection not in range, 1 if it intersects at a vertex, 2 if it intersects regularly, 3 if error(never actually happens)
{
    if (l1.m != l2.m) {
        double x_intersect =
                l2.m == INF || l1.m == INF ? (l2.m == INF ? 0 - l2.c : 0 - l1.c) : r((l2.c - l1.c) / (l1.m - l2.m));
        double y_intersect = l2.m == INF ? r(x_intersect * l1.m + l1.c) : r(x_intersect * l2.m + l2.c);

        lines_intersection[0] = x_intersect;
        lines_intersection[1] = y_intersect;

        //Checks if the line intersects the edge within the 2 end vertices
        if (!((l2.x1 == l2.x2 || l2.y1 == l2.y2 ? (l2.x1 == l2.x2 ? ((y_intersect >= (l2.y1 > l2.y2 ? l2.y2 : l2.y1)) &&
                                                                    (y_intersect <= (l2.y1 > l2.y2 ? l2.y1 : l2.y2))) : (
                                                          (x_intersect >= (l2.x1 > l2.x2 ? l2.x2 : l2.x1)) &&
                                                          (x_intersect <= (l2.x1 > l2.x2 ? l2.x1 : l2.x2)))) : (
                      x_intersect >= (l2.x1 > l2.x2 ? l2.x2 : l2.x1) && x_intersect <= (l2.x1 > l2.x2 ? l2.x1 : l2.x2))) && (l1.x1 == l1.x2 || l1.y1 == l1.y2 ? (l1.x1 == l1.x2 ? ((y_intersect >= (l1.y1 > l1.y2 ? l1.y2 : l1.y1)) &&
                                                                                                                                                                                   (y_intersect <= (l1.y1 > l1.y2 ? l1.y1 : l1.y2))) : (
                                                                                                                                                                                       (x_intersect >= (l1.x1 > l1.x2 ? l1.x2 : l1.x1)) &&
                                                                                                                                                                                       (x_intersect <= (l1.x1 > l1.x2 ? l1.x1 : l1.x2)))) : (
                                                                                                                                                   x_intersect >= (l1.x1 > l1.x2 ? l1.x2 : l1.x1) && x_intersect <= (l1.x1 > l1.x2 ? l1.x1 : l1.x2)))))
            return 0;

        for (int i = 0; i < num_vertices; ++i)
            if (x_intersect == vertices[i][0] && y_intersect == vertices[i][1]) {
                intersection_vertex_number = i;
                return 1;
            }
        return 2;
    }
    return 3;
}



_Bool point_in_polygon(double x, double y)
{
    struct Line_Double current_line = {.m = 0, .c = y, .x1 = x, .y1 = y, .x2 = INF, .y2 = y};

    int result;
    int counter = 0;
    int vertex_intersection_counter = 0;
    double val1, val2;

    for (int i = 0; i < num_edges; ++i) {

        if (current_line.m == edge_equations[i].m) {
            if (current_line.c == edge_equations[i].c) {
                int neighbours[2][2];
                int neighbour_counter = 0;
                for (int j = 0; j < num_edges; ++j) {       //Find neighbours of the vertices on the end points of the current line
                    if (j == i)
                        continue;
                    if ((edge_equations[j].x1 == edge_equations[i].x1) && (edge_equations[j].y1 == edge_equations[i].y1)) {
                        neighbours[neighbour_counter][0] = edge_equations[j].x2;
                        neighbours[neighbour_counter++][1] = edge_equations[j].y2;
                    }
                    else if ((edge_equations[j].x2 == edge_equations[i].x1) && (edge_equations[j].y2 == edge_equations[i].y1)) {
                        neighbours[neighbour_counter][0] = edge_equations[j].x1;
                        neighbours[neighbour_counter++][1] = edge_equations[j].y1;
                    }
                    else if ((edge_equations[j].x1 == edge_equations[i].x2) && (edge_equations[j].y1 == edge_equations[i].y2)) {
                        neighbours[neighbour_counter][0] = edge_equations[j].x2;
                        neighbours[neighbour_counter++][1] = edge_equations[j].y2;
                    }
                    else if ((edge_equations[j].x2 == edge_equations[i].x2) && (edge_equations[j].y2 == edge_equations[i].y2)) {
                        neighbours[neighbour_counter][0] = edge_equations[j].x1;
                        neighbours[neighbour_counter++][1] = edge_equations[j].y1;
                    }
                }

                val1 = current_line.m == INF ? (current_line.x1 - neighbours[0][0]) : (neighbours[0][1] -
                                                                                       (current_line.m *
                                                                                        neighbours[0][0]) -
                                                                                       current_line.c);
                val2 = current_line.m == INF ? (current_line.x1 - neighbours[1][0]) : (neighbours[1][1] -
                                                                                       (current_line.m *
                                                                                        neighbours[1][0]) -
                                                                                       current_line.c);

                if (!((val1 > 0 && val2 > 0) || (val1 < 0 && val2 < 0)))
                    ++counter;
            }

        } else {
            result = find_line_intersection(current_line, edge_equations[i]);

#ifdef Test
            if (result)
                printf("");
#endif

            switch (result) {
                case 2: {
                    counter++;
                    break;
                }
                case 1: {


                    //Check if both neighbours of the intersection vertex are on the same side of the line, ignore it.
                    int neighbours[2][2];
                    int neighbour_count = 0;
                    for (int j = 0; j < num_edges; ++j) {
                        if (edge_equations[j].x1 == lines_intersection[0] &&
                            edge_equations[j].y1 == lines_intersection[1]) {
                            neighbours[neighbour_count][0] = edge_equations[j].x2;
                            neighbours[neighbour_count++][1] = edge_equations[j].y2;
                        } else if (edge_equations[j].x2 == lines_intersection[0] &&
                                   edge_equations[j].y2 == lines_intersection[1]) {
                            neighbours[neighbour_count][0] = edge_equations[j].x1;
                            neighbours[neighbour_count++][1] = edge_equations[j].y1;
                        }
                    }

                    val1 = current_line.m == INF ? (current_line.x1 - neighbours[0][0]) : (neighbours[0][1] -
                                                                                           (current_line.m *
                                                                                            neighbours[0][0]) -
                                                                                           current_line.c);
                    val2 = current_line.m == INF ? (current_line.x1 - neighbours[1][0]) : (neighbours[1][1] -
                                                                                           (current_line.m *
                                                                                            neighbours[1][0]) -
                                                                                           current_line.c);

                    if (!((val1 > 0 && val2 > 0) || (val1 < 0 && val2 < 0)) && val1 && val2) {
                        ++vertex_intersection_counter;
                        ++counter;
                        break;
                    }
                }
            }
        }
    }

    counter -= (vertex_intersection_counter / 2);

    return counter % 2;
}

int main(int argc, char *argv[]) {

    double longest_line = 0;

    FILE *file = fopen(argv[1], "r");

    fscanf(file, "%d", &num_vertices);
    num_edges = num_vertices;

    vertices = (int **) malloc(num_vertices * sizeof(int *));
    for (int i = 0; i < num_vertices; ++i)
        vertices[i] = (int *) malloc(2 * sizeof(int));

    edge_equations = (struct Line *) malloc(num_vertices * sizeof(struct Line));

    fscanf(file, "%d %d", *(vertices), *(vertices) + 1);   //We need to have atleast 2 vertices for an edge to be made
    for (int i = 1; i < num_vertices; ++i) {
        fscanf(file, "%d %d", *(vertices + i), *(vertices + i) + 1);

        make_edge_equation(i /*/ 2*/ - 1, vertices[(i - 1) % num_vertices][0], vertices[(i - 1) % num_vertices][1],
                           vertices[i][0], vertices[i][1]);
    }
    make_edge_equation(num_edges - 1, vertices[num_edges - 1][0], vertices[num_edges - 1][1], vertices[0][0],
                       vertices[0][1]);

#ifdef Test

    for (int m = 0; m < num_edges; ++m)
        printf("x1: %d  y1: %d | x2: %d  y2: %d | m: %lf  c: %lf\n", edge_equations[m].x1, edge_equations[m].y1,
               edge_equations[m].x2, edge_equations[m].y2, edge_equations[m].m, edge_equations[m].c);

#endif

    double x_intersect;
    double y_intersect;

    intersepts = (struct intersept *) malloc(num_edges * sizeof(struct intersept));

    int intersept_counter = 0;
    _Bool invalid;  //Set to 1 if there if an intersection is found in between the line segment

    for (int i = 0; i < num_vertices; ++i)
        for (int j = i + 1; j < num_vertices; ++j) {

#ifdef Test
            if (i == 46 && j == 50)
                printf("\n");
#endif

            intersept_counter = 0;
            invalid = 0;
            //Construct the current line we want to consider
            struct Line current_line = FindMC(vertices[i][0], vertices[i][1], vertices[j][0], vertices[j][1]);
            current_line.x1 = vertices[i][0];
            current_line.y1 = vertices[i][1];
            current_line.x2 = vertices[j][0];
            current_line.y2 = vertices[j][1];

            if ((j != i +
                      1) && !(i == 0 && j == num_vertices - 1)/*Last edge*/) //We are not looking at an edge of the polygon. If we are, we have way fewer things to check.
                if (!point_in_polygon((1.0 * ((current_line.x2 * 1 + current_line.x1 * 1000000)) / (1 + 1000000)),
                                      (1.0 * ((current_line.y2 * 1 + current_line.y1 * 1000000)) / (1 + 1000000))))
                    continue;

            for (int l = 0; l < num_edges; ++l)
                intersepts[l].x = intersepts[l].y = DBL_MAX;

            for (int k = 0; (k < num_edges); ++k)
                if (current_line.m != edge_equations[k].m) {
                    x_intersect =
                            edge_equations[k].m == INF || current_line.m == INF ? (edge_equations[k].m == INF ? 0 -
                                                                                                                edge_equations[k].c
                                                                                                              : 0 -
                                                                                                                current_line.c)
                                                                                :
                            r((edge_equations[k].c - current_line.c) / (current_line.m - edge_equations[k].m));

                    y_intersect = edge_equations[k].m == INF ? r(x_intersect * current_line.m + current_line.c) : r(
                            x_intersect * edge_equations[k].m + edge_equations[k].c);

                    //Checks if the line intersects the edge within or on the 2 end vertices
                    if (!(edge_equations[k].x1 == edge_equations[k].x2 || edge_equations[k].y1 == edge_equations[k].y2
                          ? (edge_equations[k].x1 == edge_equations[k].x2 ? ((y_intersect >= (edge_equations[k].y1 >
                                                                                              edge_equations[k].y2
                                                                                              ? edge_equations[k].y2
                                                                                              : edge_equations[k].y1)) &&
                                                                             (y_intersect <= (edge_equations[k].y1 >
                                                                                              edge_equations[k].y2
                                                                                              ? edge_equations[k].y1
                                                                                              : edge_equations[k].y2)))
                                                                          : ((x_intersect >= (edge_equations[k].x1 >
                                                                                              edge_equations[k].x2
                                                                                              ? edge_equations[k].x2
                                                                                              : edge_equations[k].x1)) &&
                                                                             (x_intersect <= (edge_equations[k].x1 >
                                                                                              edge_equations[k].x2
                                                                                              ? edge_equations[k].x1
                                                                                              : edge_equations[k].x2))))
                          : (x_intersect >= (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x2
                                                                                         : edge_equations[k].x1) &&
                             x_intersect <= (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x1
                                                                                         : edge_equations[k].x2)/* && y_intersect > (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y2 : edge_equations[k].y1) && y_intersect < (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y1 : edge_equations[k].y2)*/)))
                        continue;

                    //Checks if the line drawn between 2 vertices intersects any other vertex or edge
                    if (current_line.m == INF ? (
                            y_intersect > (vertices[i][1] > vertices[j][1] ? vertices[j][1] : vertices[i][1]) &&
                            y_intersect < (vertices[i][1] > vertices[j][1] ? vertices[i][1] : vertices[j][1])) : (
                                x_intersect > (vertices[i][0] > vertices[j][0] ? vertices[j][0] : vertices[i][0]) &&
                                x_intersect < (vertices[i][0] > vertices[j][0] ? vertices[i][0] : vertices[j][0]))) {
                        invalid = 1;
                        break;
                    }

                    intersepts[intersept_counter].x = x_intersect;
                    intersepts[intersept_counter++].y = y_intersect;
                }

            if (!invalid) {

                double current_line_length_possible = find_length_using_intersepts(intersept_counter, vertices[i][0],
                                                                                   vertices[i][1], vertices[j][0],
                                                                                   vertices[j][1], current_line);


                longest_line =
                        current_line_length_possible > longest_line ? current_line_length_possible : longest_line;
            }
        }

    printf("Longest line in the polygon is %lf\n", longest_line);
    return 0;
}
