#ifndef ellips_h
#define ellips_h

typedef struct {
  float minor;  /* Minor axis diameter. */
  float major;  /* Major axis diameter. */
  float pa;     /* Major axis position angle measured in a clockwise sense */
                /* wrt the +ve Y axis (radians) */
  float xc;     /* The X-axis position of the center of the ellipse */
  float yc;     /* The Y-axis position of the center of the ellipse */
  float xwid;   /* X-axis extent of rectangular area enclosing ellipse. */
  float ywid;   /* Y-axis extent of rectangular area enclosing ellipse. */  
} Ellipse;

void el_define(Ellipse *el, float minor, float major, float pa,
	       float xc, float yc);

void el_move(Ellipse *el, float xc, float yc);

void el_locus(Ellipse *el, float theta, float *x, float *y);

typedef enum {EL_FULL, EL_PART, EL_CENT} Elstat;

int el_visible(Ellipse *el, float xa, float xb, float ya, float yb,
	       Elstat state);

void el_plot(Ellipse *el, int outline, int fill, int cross, int nmax);

#endif
