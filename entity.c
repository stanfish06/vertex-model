#include "entity.h"

Vertex* create_vertex(double x, double y) {
    Vertex* v = malloc(sizeof(Vertex));
    if (!v) return NULL;
    v->x = x;
    v->y = y;
    v->out = NULL;
    return v;
}

Edge* create_edge(Vertex* v) {
    Edge* e = malloc(sizeof(Edge));
    e->root = v;
    v->out = e;
    e->next = e->prev = e->cell = e->twin = NULL;
    return e;
}

Cell* create_cell(int n_faces, double* xx, double* yy) {
    Cell* c = malloc(sizeof(Cell));
    Vertex* v0 = create_vertex(xx[0], yy[0]);
    Edge* e0 = create_edge(v0);
    c->face = e0;
    e0->cell = c;
    Edge* e_prev = e0;
    Edge* e;
    for (int i = 1; i < n_faces; i++) {
	Vertex* v = create_vertex(xx[i], yy[i]);
	e = create_edge(v);
	e_prev->next = e;
	e->prev = e_prev;
	e->cell = c;
	e_prev = e;
    }
    e->next = e0;
    e0->prev = e;

    return c;
}

Cell* join_cells(Cell* c0, Cell* c1, int i_c0_e, int i_c1_e) {
    Edge* e_join_c0 = c0->face;
    Edge* e_join_c1 = c1->face;
    for (int i = 0; i < i_c0_e; i++) {
	e_join_c0 = e_join_c0->next;
    }
    for (int i = 0; i < i_c1_e; i++) {
	e_join_c1 = e_join_c1->next;
    }
    e_join_c0.twin = e_join_c1;
    e_join_c1.twin = e_join_c0;
}

Cell* create_lattice(int n_cells) {
    // TODO: voronoi diagram initialization
}
