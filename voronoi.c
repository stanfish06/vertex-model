#include "voronoi.h"
#include <stdlib.h>

Voronoi *create_voronoi(int n, double *xx, double *yy) {
  Voronoi *v = malloc(sizeof(Voronoi));
  v->num_sites = n;
  Site *sites = malloc(sizeof(Site) * n);
  for (int i = 0; i < n; i++) {
    sites[i].x = xx[i];
    sites[i].y = yy[i];
  }
  v->sites = sites;
  v->root = NULL;
	EventQueue *eq = malloc(sizeof(EventQueue));
	// sweep line will move from left to right
	// sites with minimal x-coordinates will be added first (ordered by y)

  return v;
}

int make_event(Event *e, double x, int type, Site *site, Arc *arc, double vertex_x, double vertex_y) {
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
		} else if (e_x == ptr_x) {
			// break tie with y
			if (e_type == 0) {
				double e_y = e->data.site_event.site->y;
				if (ptr_type == 0) {
					double ptr_y = ptr->data.site_event.site->y;
					if (e_y < ptr_y) {
						e->next = ptr;
					}
				}
			}
		}
		ptr_prev = ptr;
		ptr = ptr->next;
	}
	if (ptr_prev != NULL) {
		ptr_prev->next = e;
	}
	return 0;
}

// Node *add_node(bool is_arc, Node *left, Node *right) {
//   if (is_arc) {
//     Arc *data = malloc(sizeof(Arc));
//   } else {
//     Breakpoint *data = malloc(sizeof(Breakpoint));
//   }
//   Node *n = malloc(sizeof(Node)) n->is_arc = is_arc;
//   n->data = data;
//   n->left = left;
//   n->right = right;
//
//   left->right = n;
//   right->left = n;
//   return n;
// }
//
// Cell *create_lattice(int n_cells) {
//   // TODO: voronoi diagram initialization
//   double width = 10;
//   double height = 5;
// }
