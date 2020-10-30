#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define INF INT32_MAX

/*
 * Check strict math settings
 * We rely on the assumption that we store the vertices in an order where there is an edge from vertex i to i + 1
 */

_Bool point_in_polygon(int x, int y);

struct Line {
    double m;
    double c;
    int x1;
    int y1;
    int x2;
    int y2;
}*edge_equations;

struct intersept{
    double x;
    double y;
} *intersepts;

int **vertices;
int num_vertices;
int num_edges;

struct Line FindMC(int x1, int y1, int x2, int y2) {

    struct Line line = { .m = (1.0 * (y2 - y1)) / (x2 - x1), .c = y1 - (x1 * (1.0 * (y2 - y1) / 1.0 * (x2 - x1)))};
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

int is_vertex(double x, double y) {

    for (int i = 0; i < num_vertices; ++i)
        if((double)vertices[i][0] == x && (double)vertices[i][1] == y)
            return i;

    return -1;
}

_Bool consider_point(struct intersept in, int x1, int y1, int x2, int y2, double slope) {

    int index;
    if ((index = is_vertex(in.x, in.y)) > -1) {
        if (in.x == x1 && in.y == y1) {
            if (!point_in_polygon(((x2 * 1) - (x1 * 1000000)) / (1000000 - 1),
                                  ((y2 * 1) - (y1 * 1000000)) / (1000000 - 1))) {
                return 1;
            } else
                return 0;
        }
        if (in.x == x2 && in.y == y2) {
            if (!point_in_polygon(((x1 * 1) - (x2 * 1000000)) / (1000000 - 1),
                                  ((y1 * 1) - (y2 * 1000000)) / (1000000 - 1))) {
                return 1;
            } else
                return 0;
        }

        if (edge_equations[index].m == slope || edge_equations[index + 1].m == slope || edge_equations[index - 1].m == slope) {
            if (edge_equations[index].m == slope && (sqrt(pow(x2 - in.x, 2) + pow(y2 - in.y, 2)) < sqrt(pow(x2 - (edge_equations[index].x1 == in.x && edge_equations[index].y1 == in.y ? edge_equations[index].x2 : edge_equations[index].x1), 2) + pow(y2 - (edge_equations[index].x1 == in.x && edge_equations[index].y1 == in.y ? edge_equations[index].y2 : edge_equations[index].y1), 2)))) {
                return 0;
            }
            else if (edge_equations[index + 1].m == slope && (sqrt(pow(x2 - in.x, 2) + pow(y2 - in.y, 2)) < sqrt(pow(x2 - (edge_equations[index + 1].x1 == in.x && edge_equations[index + 1].y1 == in.y ? edge_equations[index + 1].x2 : edge_equations[index + 1].x1), 2) + pow(y2 - (edge_equations[index + 1].x1 == in.x && edge_equations[index + 1].y1 == in.y ? edge_equations[index + 1].y2 : edge_equations[index + 1].y1), 2)))) {
                return 0;
            }
            else if (edge_equations[index - 1].m == slope && (sqrt(pow(x2 - in.x, 2) + pow(y2 - in.y, 2)) < sqrt(pow(x2 - (edge_equations[index - 1].x1 == in.x && edge_equations[index - 1].y1 == in.y ? edge_equations[index - 1].x2 : edge_equations[index - 1].x1), 2) + pow(y2 - (edge_equations[index - 1].x1 == in.x && edge_equations[index - 1].y1 == in.y ? edge_equations[index - 1].y2 : edge_equations[index - 1].y1), 2)))) {
                return 0;
            }
            int point1[2], point2[2];
            for (int i = index - 1; i < index + 2; ++i) {
                if (edge_equations[i].m != slope) {
                    if(edge_equations[i].x1 == in.x && edge_equations[i].y1 == in.y) {
                        point1[0] = edge_equations[i].x2;
                        point1[1] = edge_equations[i].y2;
                    }
                }
            }
        }


    }

    return 1;
}


double find_length_using_intersepts(int intersept_counter, int x1, int y1, int x2, int y2, double slope) {

    double sideA = DBL_MAX, sideB = DBL_MAX;

    for (int i = 0; i < intersept_counter; ++i) {

        if (consider_point(intersepts[i], x1, y1, x2, y2, slope)) {
                double d1 = sqrt(pow((intersepts[i].y - y1), 2) + pow((intersepts[i].x - x1), 2));
                double d2 = sqrt(pow((intersepts[i].y - y2), 2) + pow((intersepts[i].x - x2), 2));

                if (d1 == 0.0)
                    sideA = 0;
                if (d2 == 0.0)
                    sideB = 0;

                if (d1 > d2) {
                    if (sideA > d1)
                        sideA = d1;
                }
                else if (d2 > d1) {
                    if (sideB > d2)
                        sideB = d2;
                }
        }
    }

    //sideA = sideA == DBL_MAX ? 0 : sideA;
    //sideB = sideB == DBL_MAX ? 0 : sideB;
}


int find_line_intersection(struct Line l1, struct Line l2) //0 if intersection not in range, 1 if it intersects at a vertex, 2 if it intersects regularly, 3 if error(never actually happens)
{
    if (l1.m != l2.m) {
        double x_intersect = l2.m == INF || l1.m == INF ? (l2.m == INF ? 0 - l1.c : 0 - l2.c) : (l2.c - l1.c) / (l1.m - l2.m);
        double y_intersect = l2.m == INF ? (x_intersect * l1.m + l1.c) : (x_intersect * l2.m + l2.c);

        lines_intersection[0] = x_intersect;
        lines_intersection[1] = y_intersect;

        //Checks if the line intersects the edge within the 2 end vertices
        if (!(l2.x1 == l2.x2 || l2.y1 == l2.y2 ? (l2.x1 == l2.x2 ? ((y_intersect > (l2.y1 > l2.y2 ? l2.y2 : l2.y1)) && (y_intersect < (l2.y1 > l2.y2 ? l2.y1 : l2.y2))) : ((x_intersect > (l2.x1 > l2.x2 ? l2.x2 : l2.x1)) && (x_intersect < (l2.x1 > l2.x2 ? l2.x1 : l2.x2)))) : (x_intersect > (l2.x1 > l2.x2 ? l2.x2 : l2.x1) && x_intersect < (l2.x1 > l2.x2 ? l2.x1 : l2.x2))))
            return 0;

        for (int i = 0; i < num_vertices; ++i)
            if(x_intersect == vertices[i][0] && y_intersect == vertices[i][1]) {
                intersection_vertex_number = i;
                return 1;
            }
        return 2;
    }
    return 3;
}



_Bool point_in_polygon(int x, int y)
{
    struct Line current_line = {.m = 0, .c = y, .x1 = x, .y1 = y, .x2 = INT32_MAX, .y2 = y};

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

            switch (result) {
                case 1: {
                    counter++;
                    break;
                }
                case 2: {

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

                    if (!((val1 > 0 && val2 > 0) || (val1 < 0 && val2 < 0))) {
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

int main() {

    double longest_line = 0;

    scanf("%d", &num_vertices);
    num_edges = num_vertices;

    vertices = (int **) malloc(num_vertices * sizeof(int *));
    for (int i = 0; i < num_vertices; ++i)
        vertices[i] = (int *) malloc(2 * sizeof(int));

    edge_equations = (struct Line *) malloc(num_vertices * sizeof(struct Line));

    for (int i = 0; i < num_vertices; ++i) {
        scanf("%d %d", *(vertices + i), *(vertices + i) + 1);

        //if (i % 2 && i)     //TODO: Fix later
        make_edge_equation(i / 2 - 1, vertices[i - 1][0], vertices[i - 1][1], vertices[i][0], vertices[i][1]);
    }

    double x_intersect;
    double y_intersect;

    intersepts = (struct intersept *) malloc(num_edges * sizeof(struct intersept));

    int intersept_counter = 0;
    _Bool invalid;  //Set to 1 if there if an intersection is found in between the line segment

    for (int i = 0; i < num_vertices; ++i)
        for (int j = i + 1; j < num_vertices; ++j) {

            invalid = 0;
            //Construct the current line we want to consider
            struct Line current_line = FindMC(vertices[i][0], vertices[i][1], vertices[j][0], vertices[j][1]);

            if (j != i + 1) //We are not looking at an edge of the polygon. If we are, we have way fewer things to check.
                if (!point_in_polygon(((current_line.x2 * 1 + current_line.x1 * 1000000) / (1 + 1000000)), ((current_line.y2 * 1 + current_line.y1 * 1000000) / (1 + 1000000))))
                    continue;


            for (int l = 0; l < num_edges; ++l)
                intersepts[l].x = intersepts[l].y = DBL_MAX;

            for (int k = 0; (k < num_edges); ++k)
                if (current_line.m != edge_equations[k].m) {
                    x_intersect = edge_equations[k].m == INF || current_line.m == INF ? (edge_equations[k].m == INF ? 0 - current_line.c : 0 - edge_equations[k].c) : (edge_equations[k].c - current_line.c) / (current_line.m - edge_equations[k].m);
                    y_intersect = edge_equations[k].m == INF ? (x_intersect * current_line.m + current_line.c) : (x_intersect * edge_equations[k].m + edge_equations[k].c);

                    //Checks if the line intersects the edge within the 2 end vertices (possibly wrong)
                    //if (!(edge_equations[k].x1 == edge_equations[k].x2 || edge_equations[k].y1 == edge_equations[k].y2 ? (edge_equations[k].x1 == edge_equations[k].x2 ? ((y_intersect > (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y2 : edge_equations[k].y1)) && (y_intersect < (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y1 : edge_equations[k].y2))) : ((x_intersect > (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x2 : edge_equations[k].x1)) && (x_intersect < (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x1 : edge_equations[k].x2)))) : (x_intersect > (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x2 : edge_equations[k].x1) && x_intersect < (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x1 : edge_equations[k].x2)/* && y_intersect > (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y2 : edge_equations[k].y1) && y_intersect < (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y1 : edge_equations[k].y2)*/)))
                    //Checks if the line intersects the edge within or on the 2 end vertices
                    if (!(edge_equations[k].x1 == edge_equations[k].x2 || edge_equations[k].y1 == edge_equations[k].y2 ? (edge_equations[k].x1 == edge_equations[k].x2 ? ((y_intersect >= (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y2 : edge_equations[k].y1)) && (y_intersect <= (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y1 : edge_equations[k].y2))) : ((x_intersect >= (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x2 : edge_equations[k].x1)) && (x_intersect <= (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x1 : edge_equations[k].x2)))) : (x_intersect >= (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x2 : edge_equations[k].x1) && x_intersect <= (edge_equations[k].x1 > edge_equations[k].x2 ? edge_equations[k].x1 : edge_equations[k].x2)/* && y_intersect > (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y2 : edge_equations[k].y1) && y_intersect < (edge_equations[k].y1 > edge_equations[k].y2 ? edge_equations[k].y1 : edge_equations[k].y2)*/)))
                        continue;

                    //Checks if the line drawn between 2 vertices intersects any other vertex or edge
                    if (current_line.m == INF ? (y_intersect >= (vertices[i][1] > vertices[j][1] ? vertices[j][1] : vertices[i][1]) && y_intersect <= (vertices[i][1] > vertices[j][1] ? vertices[i][1] : vertices[j][1])) : (x_intersect >= (vertices[i][0] > vertices[j][0] ? vertices[j][0] : vertices[i][0]) && x_intersect <= (vertices[i][0] > vertices[j][0] ? vertices[i][0] : vertices[j][0]))) {
                        invalid = 1;
                        break;
                    }

                    intersepts[intersept_counter].x = x_intersect;
                    intersepts[intersept_counter++].y = y_intersect;
                }

            if (!invalid) {
                double current_line_length_possible = find_length_using_intersepts(intersept_counter, vertices[i][0], vertices[i][1], vertices[j][0], vertices[j][1], current_line.m);
                longest_line = current_line_length_possible > longest_line ? current_line_length_possible : longest_line;
            }
        }

    printf("Longest line in the polygon is %lf\n", longest_line);
    return 0;
}