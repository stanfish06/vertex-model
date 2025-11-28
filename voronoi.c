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

Circle *compute_circumcircle(Site *p1, Site *p2, Site *p3) {
  Circle *c = malloc(sizeof(Circle));
  double x1 = p1->x;
  double y1 = p1->y;
  double x2 = p2->x;
  double y2 = p2->y;
  double x3 = p3->x;
  double y3 = p3->y;

  // midpoint between p1 and p2
  double m1x = (x1 + x2) / 2;
  double m1y = (y1 + y2) / 2;
  // midpoint between p2 and p3
  double m2x = (x2 + x3) / 2;
  double m2y = (y2 + y3) / 2;

  // slope between p1 and p2
  double s1 = (y2 - y1) / (x2 - x1);
  // slope between p2 and p3
  double s2 = (y3 - y2) / (x3 - x2);

  // solve m1y - 1/s1(x-m1x) = m2y - 1/s2(x-m2x) -> x = (m2y - m1y + m2x/s2 -
  // m1x/s1)/(1/s2 - 1/s1)
  double cx = (m2y - m1y + m2x / s2 - m1x / s1) / (1 / s2 - 1 / s1);
  double cy = m1y - 1 / s1 * (cx - m1x);
  double r = sqrt(pow(x1 - cx, 2) + pow(y1 - cy, 2));
  c->cx = cx;
  c->cy = cy;
  c->r = r;
  return c;
}

/*
 * after adding a new site, check two incident triplets to see if breakpoints
will merge:
 * 1. for each triplet, find circumcircle
 * 2. since my sweep line sweeps from left to right, check if the right most
point of the circle is on the left side of the sweep line
 * 3. if so, add circle event to queue
 */
void check_circle_after_site_event(Node *n) {
  Node *ptr = n;
  Site *p1 = (*n).data.arc->site;
  Site *p2 = NULL;
  Site *p3 = NULL;
  // up triplet
  ptr = ptr->left;
  while ((p2 == NULL || p3 == NULL) && ptr != NULL) {
    if (ptr->is_arc) {
      if (p2 == NULL) {
        p2 = (*ptr).data.arc->site;
      } else {
        p3 = (*ptr).data.arc->site;
      }
    }
    ptr = ptr->left;
  }
  Circle *c = compute_circumcircle(p1, p2, p3);
  // current position of sweep line is just the newly added node
  if (c->cx + c->r > p1->x) {
    // add circle event
  }
  free(c);

  // down triplet
  p2 = NULL;
  p3 = NULL;
  ptr = ptr->right;
  while ((p2 == NULL || p3 == NULL) && ptr != NULL) {
    if (ptr->is_arc) {
      if (p2 == NULL) {
        p2 = (*ptr).data.arc->site;
      } else {
        p3 = (*ptr).data.arc->site;
      }
    }
    ptr = ptr->right;
  }
  c = compute_circumcircle(p1, p2, p3);
  // current position of sweep line is just the newly added node
  if (c->cx + c->r > p1->x) {
    // add circle event
  }
}

/*
 * implement this as well
 * the position of sweep line is just the right most x of the circumcircle as it
 * is tangent to that
 */
void check_circle_after_circle_event(Node *n) {}

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
    // TODO: need to free data as well
    free(e);
  }
  free(eq);
  return v;
}
