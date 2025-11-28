#include "voronoi.h"
#include <stdio.h>

void testing_eventqueue(int n, double *xx, double *yy) {
  EventQueue *eq = malloc(sizeof(EventQueue));
  eq->first = NULL;
  for (int i = 0; i < n; ++i) {
    Site *s = malloc(sizeof(Site));
    Event *e = malloc(sizeof(Event));
    s->x = xx[i];
    s->y = yy[i];
    make_event(e, s->x, 0, s, 0, 0, 0);
    add_event(eq, e);
  }
  Event *e;
  while (eq->first) {
    e = pop_event(eq);
    char *type = (e->type == 0) ? "site event" : "circle event";
    double x = e->x;
    printf("%s\n", type);
    printf("x=%f\n", x);
    Site *s;
    if (e->type == 0) {
      s = (*e).data.site_event.site;
      printf("(x,y)=(%f,%f)\n", s->x, s->y);
      free(s);
    } else {
      printf("Circle event\n");
    }
    free(e);
  }
  free(eq);
}

int main() {
  double xx[] = {1.0, 2.0, 2.5, 2.5, 3.0, 3.25};
  double yy[] = {0.0, -0.5, 2.0, 1.0, 5.0, 0.5};
  testing_eventqueue(6, xx, yy);
  return 0;
}
