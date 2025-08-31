#ifndef ENTITY_H
#define ENTITY_H
typedef struct Vertex Vertex; 
typedef struct Edge Edge; 
typedef struct Cell Cell; 

struct Vertex {
    double x, y;
    Edge *out;
};

struct Edge {
    Vertex *root;
    Edge *twin;
    Edge *next;
    Edge *prev;
    Cell *cell;
    double length;
};

struct Cell {
    Edge *face;
    double area;
    double perimeter;
};
#endif
