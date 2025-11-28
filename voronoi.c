#include "voronoi.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

void print_inorder(struct Node *node) {
  if (!node)
    return;
  print_inorder(node->left);
  if (node->is_arc)
    printf("Arc of the site at (%f, %f)\n", node->data.arc->site->x,
           node->data.arc->site->y);
  if (!node->is_arc)
    printf("Breakpoint between up (%f, %f) and down (%f, %f) \n",
           node->data.breakpoint->up->x, node->data.breakpoint->up->y,
           node->data.breakpoint->down->x, node->data.breakpoint->down->y);
  print_inorder(node->right);
}

int make_event(Event *e, double x, int type, Site *site, Arc *arc,
               double vertex_x, double vertex_y) {
  struct Event etmp = {0};
  etmp.x = x;
  if (type == 0) {
    etmp.data.site_event.site = site;
  } else if (type == 1) {
    etmp.data.circle_event.arc = arc;
    etmp.data.circle_event.vertex_x = vertex_x;
    etmp.data.circle_event.vertex_y = vertex_y;
  }
  etmp.type = type;
  *e = etmp;
  return 0;
}

int add_event(EventQueue *eq, Event *e) {
  Event *ptr = eq->first;
  Event *ptr_prev = NULL;
  int e_type = e->type;
  double e_x = e->x;
  while (ptr != NULL) {
    double ptr_x = ptr->x;
    int ptr_type = ptr->type;
    if (e_x < ptr_x) {
      e->next = ptr;
      break;
    } else if (e_x == ptr_x) {
      // break tie with y
      if (e_type == 0) {
        double e_y = e->data.site_event.site->y;
        if (ptr_type == 0) {
          double ptr_y = ptr->data.site_event.site->y;
          if (e_y < ptr_y) {
            e->next = ptr;
            break;
          }
        }
      }
    }
    ptr_prev = ptr;
    ptr = ptr->next;
  }
  if (ptr_prev != NULL) {
    ptr_prev->next = e;
  } else {
    // if event e is the new root
    eq->first = e;
  }
  return 0;
}

Event *pop_event(EventQueue *eq) {
  Event *e = eq->first;
  eq->first = eq->first->next;
  return e;
}

Node *add_node(bool is_arc, Node *left, Node *right) {
  Node *n = malloc(sizeof(Node));
  if (is_arc) {
    Arc *data = malloc(sizeof(Arc));
    Site *s = malloc(sizeof(Site));
    n->data.arc = data;
    n->data.arc->site = s;
  } else {
    Site *s_up = malloc(sizeof(Site));
    Site *s_down = malloc(sizeof(Site));
    Breakpoint *data = malloc(sizeof(Breakpoint));
    n->data.breakpoint = data;
    n->data.breakpoint->up = s_up;
    n->data.breakpoint->down = s_down;
  }
  n->is_arc = is_arc;
  n->left = left;
  n->right = right;

  return n;
}

// mgight be incorrect, check later
double calculate_breakpoint_x(Site *site1, Site *site2, double sweep_x) {
  double x1 = site1->x, y1 = site1->y;
  double x2 = site2->x, y2 = site2->y;

  if (fabs(y1 - y2) < 1e-9) {
    return (x1 + x2) / 2.0;
  }

  double denom = y1 - y2;
  double a = 1.0 / denom;
  double b = -2.0 * (x1 * (y2 - sweep_x) - x2 * (y1 - sweep_x)) / denom;
  double c =
      ((x1 * x1 - x2 * x2) * (y1 - y2) + (y1 * y1 - y2 * y2) * (sweep_x - y1) +
       (y2 * y2 - y1 * y1) * (sweep_x - y2)) /
      denom;

  double discriminant = b * b - 4 * a * c;
  if (discriminant < 0)
    return NAN;

  double x_plus = (-b + sqrt(discriminant)) / (2 * a);
  double x_minus = (-b - sqrt(discriminant)) / (2 * a);

  return (y1 < y2) ? x_minus : x_plus;
}

typedef struct {
  bool is_left;
  Node *parent;
  Node *child;
} Edge;
// traverse down the tree until finding an arc
// check y relative to breakpoint and choose which subtree to goto
Edge *find_arc_left(Voronoi *v, double y, double sweep_line_x) {
  Node *ptr = v->root;
  Node *ptr_prev = NULL;
  bool is_left = false;
  while (!ptr->is_arc) {
    Breakpoint *p = ptr->data.breakpoint;
    double px = calculate_breakpoint_x(p->up, p->down, sweep_line_x);
    double p1x = p->up->x;
    double p1y = p->up->y;
    double p2x = p->down->x;
    double p2y = p->down->y;
    // first you need to calculate current x based on sweep_x
    // y = 0.5*[(2x - x1 - x2)(x1 - x2)/(y2 - y1) + (y1 + y2)]
    double py =
        0.5 * ((2 * px - p1x - p2x) * (p1x - p2x) / (p2y - p1y) + (p1y + p2y));

    ptr_prev = ptr;
    if (y >= py) {
      ptr = ptr->left;
      is_left = true;
    } else {
      ptr = ptr->right;
      is_left = false;
    }
  }
  Edge *edge = malloc(sizeof(Edge));
  edge->is_left = is_left;
  edge->parent = ptr_prev;
  edge->child = ptr;
  return edge;
}

// steps:
// 1. do in-order traversal
//    - no need for full traversal (not sure what do check for now)
// 2. for each three adjacent arcs, compute their circumcircle
// 3. if exists, then find right most position of the cirlce
// 4. add circle event to queue for that position
//    - check pror circle events invovling the arc that will be deleted
// 5. process that event:
//    - remove middle arc
//    - merge two breakpoints
//    - store that voronoi vertex
// 6. recheck circle events after processing
Event *check_cycle(Voronoi *v) { return NULL; }

void process_event(Voronoi *v, Event *e) {
  int type = e->type;
  Site *site = e->data.site_event.site;
  if (type == 0) {
    if (v->root == NULL) {
      Arc *arc = malloc(sizeof(Arc));
      arc->site = site;
      Node *node = add_node(true, NULL, NULL);
      node->data.arc = arc;
      v->root = node;
      return;
    }
    double x = e->data.site_event.site->x;
    double y = e->data.site_event.site->y;

    Edge *edge = find_arc_left(v, y, x);
    Node *breakpoint_node = edge->parent;
    bool is_left = edge->is_left;
    Node *arc_node = edge->child;
    double sx = arc_node->data.arc->site->x;
    double sy = arc_node->data.arc->site->y;
    // knowing y is the new site to add,
    // this statement computes the x position where the new arc (horizontal line
    // for now) intersects with the left arc. not super sure if this is going to
    // work.
    double px = calculate_breakpoint_x(arc_node->data.arc->site, site, x);
    Breakpoint *p1 = malloc(sizeof(Breakpoint));
    Breakpoint *p2 = malloc(sizeof(Breakpoint));
    double p12x = px;
    double p12y = pow(p12x - sx, 2) / (2 * (sy - x)) + (sy + x) / 2;
    p1->x = p12x;
    p2->x = p12x;
    Arc *up = malloc(sizeof(Arc));
    Arc *middle = malloc(sizeof(Arc));
    Arc *down = malloc(sizeof(Arc));
    p1->up = arc_node->data.arc->site;
    p1->down = e->data.site_event.site;
    p2->up = e->data.site_event.site;
    p2->down = arc_node->data.arc->site;

    up->site = arc_node->data.arc->site;
    middle->site = e->data.site_event.site;
    down->site = arc_node->data.arc->site;
    // update tree
    Node *up_node = add_node(true, NULL, NULL);
    Node *middle_node = add_node(true, NULL, NULL);
    Node *down_node = add_node(true, NULL, NULL);
    // add three new arcs
    up_node->data.arc = up;
    middle_node->data.arc = middle;
    down_node->data.arc = down;
    // add two new break points
    //      p2_node (breakpoint)
    //     /                    \
    // p1_node (breakpoint)    down_node (arc)
    // /                \
    // up_node (arc)    middle_node (arc)
    Node *p1_node;
    Node *p2_node;
    p1_node = add_node(false, up_node, middle_node);
    p2_node = add_node(false, p1_node, down_node);
    if (breakpoint_node == NULL) {
      v->root = p2_node; // arc was root
    } else if (is_left) {
      breakpoint_node->left = p2_node;
    } else {
      breakpoint_node->right = p2_node;
    }
    p1_node->data.breakpoint = p1;
    p2_node->data.breakpoint = p2;
  }
}

/* ------------------------------
 * sweep from left to right
 * - site events should be added first
 * - circle events added in the main loop
 * ------------------------------ */
Voronoi *generate_voronoi(int n, double *xx, double *yy) {
  Site *sites = malloc(sizeof(Site) * n);
  EventQueue *eq = malloc(sizeof(EventQueue));
  Voronoi *v = malloc(sizeof(Voronoi));

  // sites with minimal x-coordinates will be added first (ordered by y)
  for (int i = 0; i < n; i++) {
    double x = xx[i];
    double y = yy[i];

    sites[i].x = x;
    sites[i].y = y;

    Event *e = malloc(sizeof(Event));
    make_event(e, x, 0, &sites[i], NULL, 0, 0);
    add_event(eq, e);
  }

  v->num_sites = n;
  v->sites = sites;
  v->root = NULL;

  // main loop
  while (eq->first != NULL) {
    Event *e = pop_event(eq);
    process_event(v, e);
    free(e);
  }
  free(eq);
  return v;
}
