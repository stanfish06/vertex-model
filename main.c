#include "voronoi.h"
#include <stdio.h>

int main() {
  printf("Test:\n");
  // test
  // create voronoi
  Voronoi v;
  // create this tree
  // add two new break points
  //      p1_node (breakpoint)
  //     /                    \
  // p2_node (breakpoint)    a3
  // /                \
  // a1               a2
  Node *a1 = add_node(true, NULL, NULL);
  Node *a2 = add_node(true, NULL, NULL);
  Node *a3 = add_node(true, NULL, NULL);
  Node *p2 = add_node(false, a1, a2);
  Node *p1 = add_node(false, p2, a3);
  a1->data.arc->site->x = 1.0;
  a1->data.arc->site->y = 3.0;
  a2->data.arc->site->x = 2.0;
  a2->data.arc->site->y = 2.0;
  a3->data.arc->site->x = 1.5;
  a3->data.arc->site->y = 1.0;

  p2->data.breakpoint->up = a1->data.arc->site;
  p2->data.breakpoint->down = a2->data.arc->site;
  p1->data.breakpoint->up = a3->data.arc->site;
  p1->data.breakpoint->down = a2->data.arc->site;
  p1->data.breakpoint->x = 1.3;
  p2->data.breakpoint->x = 1.6;

  v.root = p1;
  print_inorder(v.root);

  free(a1);
  free(a2);
  free(a3);
  free(p1);
  free(p2);
  return 0;
}
