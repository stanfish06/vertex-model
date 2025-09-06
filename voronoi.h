#ifndef VORONOI_H
#define VORONOI_H
typedef struct Node Node;

typedef struct {
  double x, y;
} Site;

typedef struct {
  // left and right seed points
  // use bisector equation to get boundary ray
  Site *left;
  Site *right;
  // current position of the breakpoint
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
  Node *left;
  Node *right;
};

typedef struct {
  Node *root;
  Site *sites;
  int num_sites;
} Voronoi;

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
typedef struct {
  Event *first;
} EventQueue;
#endif
