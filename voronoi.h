#ifndef VORONOI_H
#define VORONOI_H
typedef struct Node Node;
#include <stdbool.h>
#include <stdlib.h>

typedef struct {
  double x, y;
} Site;

// To find y
// (x - x1)^2 + (y - y1)^2 = (x - x2)^2 + (y - y2)^2
// (2y - y1 - y2)(y2 - y1) = (2x - x1 - x2)(x1 - x2)
// y = 0.5*[(2x - x1 - x2)(x1 - x2)/(y2 - y1) + (y1 + y2)]
typedef struct {
  // up and down seed points
  // use bisector equation to get boundary ray
  Site *up;
  Site *down;
  // initial position of the breakpoint
  // upon disappearing, you can derive the equation for the boundary ray
  double x;
} Breakpoint;

typedef struct {
  // root of the arc
  Site *site;
} Arc;

struct Node {
  int is_arc;
  union {
    Arc *arc;
    Breakpoint *breakpoint;
  } data;
  Node *left; // left = up
  Node *right; // right = down
};

void print_inorder(struct Node *node);
Node *add_node(bool is_arc, Node *left, Node *right);

// event queue
// two events:
// add parabola
// remove arc
typedef struct Event Event;
struct Event {
  double x;
  int type;
  union {
    struct {
      Site *site;
    } site_event;
    struct {
      double vertex_x;
      double vertex_y;
      Arc *arc;
    } circle_event;
  } data;
  Event *next;
};

int make_event(Event *e, double x, int type, Site *site, Arc *arc, double vertex_x, double vertex_y);

typedef struct {
  Event *first;
} EventQueue;


typedef struct {
  Node *root;
  Site *sites;
  int num_sites;
} Voronoi;
#endif
