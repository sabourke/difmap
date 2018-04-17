#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "color.h"
#include "symtab.h"
#include "logio.h"
#include "cpgplot.h"

#define BASE 17      /* Lowest color index that PGGRAY will use. */
#define MAXLEVS 127  /* Max number of gray-scale levels in PGPLOT */
#define MINLEVS 15   /* Min number of gray-scale levels in PGPLOT */

static void *del_SCmap(void *value);
static int plcmap(float *l, float *r, float *g, float *b, int nc,
		  float contra, float bright);


/*.......................................................................
 * Define the default color tables.
 */

/*
 * Define single-color ramp functions.
 */
static float grey_l[]  = {0.0,1.0};
static float grey_c[]  = {0.0,1.0};
static float blank_c[] = {0.0,0.0};
/*
 * Define a rainbow color table.
 */
static float rain_l[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
static float rain_r[] = { 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
static float rain_g[] = { 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
static float rain_b[] = { 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};

/*
 * Iraf "heat" color table.
 */
static float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
static float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
static float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
static float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};

/*
 * Weird IRAF ramp color table.
 */
static float ramp_l[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0};
static float ramp_r[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
static float ramp_g[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
static float ramp_b[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};

/*
 * AIPS tvfiddle discrete rainbow color table.
 */
static float aips_l[] = {0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
			 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
static float aips_r[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
			 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
static float aips_g[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
			 0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
static float aips_b[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
			 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

#ifdef A_SIZE
#undef A_SIZE
#endif
#define A_SIZE(lev) sizeof(lev)/sizeof(lev[0])

static Cmap std_cmaps[] = {
  {"aips",   CM_COLOR, 1, A_SIZE(aips_l),   aips_l,  aips_r,  aips_g,  aips_b},
  {"blue"   ,CM_COLOR, 1, A_SIZE(grey_l),   grey_l, blank_c, blank_c,  grey_c},
  {"green"  ,CM_COLOR, 1, A_SIZE(grey_l),   grey_l, blank_c,  grey_c, blank_c},
  {"grey",    CM_GREY, 1, A_SIZE(grey_l),   grey_l,  grey_c,  grey_c,  grey_c},
  {"heat",   CM_COLOR, 1, A_SIZE(heat_l),   heat_l,  heat_r,  heat_g,  heat_b},
  {"none",   CM_NONE,  1, 0,                  NULL,    NULL,    NULL,    NULL},
  {"rainbow",CM_COLOR, 1, A_SIZE(rain_l),   rain_l,  rain_r,  rain_g,  rain_b},
  {"ramp",   CM_COLOR, 1, A_SIZE(ramp_l),   ramp_l,  ramp_r,  ramp_g,  ramp_b},
  {"red"    ,CM_COLOR, 1, A_SIZE(grey_l),   grey_l,  grey_c, blank_c, blank_c},
};
static int n_std_cmap = A_SIZE(std_cmaps);

/*.......................................................................
 * Describe recognised transfer functions.
 */
static struct {
  char *name;     /* The name to give the transfer function */
  Cmtran type;    /* The enumerated type of the transfer function */
} std_trans[] = {
  {"linear", TR_LINEAR},
  {"log",    TR_LOG},
  {"sqrt",   TR_SQRT}
};
static int n_std_trans = A_SIZE(std_trans);

/*.......................................................................
 * Create a color table.
 *
 * Output:
 *  return   Ctable *  The new color table, initialized with the default
 *                     color maps, or NULL on error.
 */
Ctable *new_Ctable(void)
{
  Ctable *ctab;   /* The new color table */
  int i;
/*
 * Allocate the container.
 */
  ctab = (Ctable *) malloc(sizeof(Ctable));
  if(ctab == NULL) {
    lprintf(stderr, "new_Ctable: Insufficient room for color table.\n");
    return NULL;
  };
/*
 * Initialize the container at least to the point at which it can safely be
 * sent to del_Ctable().
 */
  ctab->symtab = NULL;
  ctab->cmap = NULL;
  ctab->contra = 1.0f;
  ctab->bright = 0.5f;
  ctab->tran = TR_LINEAR;
  ctab->vmin = 0.0f;
  ctab->vmax = 0.0f;
/*
 * Allocate the new color-map symbol table with room for all the standard
 * color maps + 1 to allow us to install a changeable default color map
 * from one of the standard color maps.
 */
  ctab->symtab = new_Symtab(n_std_cmap+1, "Colormap", amb_report, del_SCmap);
  if(ctab->symtab==NULL)
    return del_Ctable(ctab);
/*
 * Insert each of the standard color maps.
 */
  for(i=0; i<n_std_cmap; i++) {
    if(add_Cmap(ctab, std_cmaps[i].name, &std_cmaps[i]))
      return i==0 ? del_Ctable(ctab) : ctab;
  };
/*
 * Install the rainbow colormap as the default "color" colormap.
 */
  ctab->cmap = get_Cmap(ctab, "rainbow");
/*
 * Make the grey-scale colormap the default colormap.
 */
  ctab->cmap = get_Cmap(ctab, "grey");
  if(ctab->cmap==NULL)
    return del_Ctable(ctab);
/*
 * Install the linear transfer function as default.
 */
  ctab->tran = get_Cmtran("linear");
/*
 * Return the initialized color table.
 */
  return ctab;
}

/*.......................................................................
 * Delete a table of color maps. Note that del_Cmap() is called on each
 * color table, to delete the symbols (only if their reference counts
 * are zero).
 *
 * Input:
 *  ctab   Ctable *   The color table to be deleted.
 * Output:
 *  return Ctable *   Allways NULL.
 */
Ctable *del_Ctable(Ctable *ctab)
{
  if(ctab) {
    ctab->symtab = del_Symtab(ctab->symtab);
    free(ctab);
  };
  return NULL;
}

/*.......................................................................
 * Define a function that can be called by del_Symtab() to delete a
 * color map.
 *
 * Input:
 *  value    void *   A colormap cast to (void *).
 * Output:
 *  return   void *   Allways NULL.
 */
static void *del_SCmap(void *value)
{
  return (void *) del_Cmap((Cmap *) value);
}

/*.......................................................................
 * Decrement the reference count of a color map and delete it if the
 * result is <= 0. This combined with the fact that add_Cmap() increments
 * the reference count, allows a single color map to appear more than
 * once in one or many color tables, and allows statically allocated
 * color maps (reference count must be initialized to 1) to be used without
 * the danger of trying to delete them.
 *
 * Input:
 *  cmap    Cmap *   The color map to be deleted.
 * Output:
 *  return  Cmap *   Allways NULL.
 */
Cmap *del_Cmap(Cmap *cmap)
{
  if(cmap && --cmap->nref <= 0) {
    if(cmap->name)
      free(cmap->name);
/*
 * All the arrays are from one array anchored at cmap->l.
 */
    if(cmap->l)
      free(cmap->l);
  };
  return NULL;
}

/*.......................................................................
 * Create a new colormap and initialize its reference count to 0.
 *
 * Input:
 *  name    char *   The name to fgive the color table. This will be
 *                   copied.
 *  nc       int     The number of colormap entries to allocate.
 * Output:
 *  return  Cmap *   The new colormap, or NULL on error.
 */
Cmap *new_Cmap(char *name, int nc)
{
  static char *nomem = "new_Cmap: Insufficient memory.\n";
  Cmap *cmap;   /* The return descriptor */
/*
 * Bad name?
 */
  if(name==NULL) {
    lprintf(stderr, "new_Cmap: NULL name intercepted.\n");
    return NULL;
  };
/*
 * Allocate the colormap container.
 */
  cmap = (Cmap *) malloc(sizeof(Cmap));
  if(cmap==NULL) {
    lprintf(stderr, nomem);
    return del_Cmap(cmap);
  };
/*
 * Initialize it at least to the point at which it is safe to send it
 * to del_Cmap().
 */
  cmap->name = NULL;
  cmap->class = CM_COLOR;
  cmap->nref = 0;
  cmap->nc = nc;
  cmap->l = NULL;
  cmap->r = NULL;
  cmap->g = NULL;
  cmap->b = NULL;
/*
 * Allocate a copy of 'name'.
 */
  cmap->name = (char *) malloc(strlen(name) + 1);
  if(cmap->name) {
    strcpy(cmap->name, name);
  } else {
    lprintf(stderr, nomem);
    return del_Cmap(cmap);
  };
/*
 * Allocate a single array of length 'nc*4' to be split between the
 * four level arrays.
 */
  cmap->l = (float *) malloc(sizeof(float) * nc * 4);
  if(cmap->l == NULL) {
    lprintf(stderr, nomem);
    return del_Cmap(cmap);
  };
/*
 * Split the array between the four arrays.
 */
  cmap->r = cmap->l + nc;
  cmap->g = cmap->r + nc;
  cmap->b = cmap->g + nc;
/*
 * Return the color map.
 */
  return cmap;
}

/*.......................................................................
 * Add a colormap to a color table.
 *
 * Input:
 *  ctab  Ctable *   The color table to add ctab to.
 *  name    char *   The name to be used to look up the color map. This
 *                   need not be the same as cmap->name.
 *  cmap    Cmap *   The color map to be added.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int add_Cmap(Ctable *ctab, char *name, Cmap *cmap)
{
  int i;
/*
 * Check that the colormap levels are arranged in increasing order.
 */
  for(i=1; i<cmap->nc; i++) {
    if(cmap->l[i] < cmap->l[i-1]) {
      lprintf(stderr, "add_Cmap: Colormap levels not in increasing order.\n");
      return 1;
    };
  };
/*
 * Add the color map to the colortable.
 */
  if(add_symbol(ctab->symtab, name, (void *) cmap, 1))
    return 1;
/*
 * Increment the reference count of the color map.
 */
  cmap->nref++;
  return 0;
}

/*.......................................................................
 * Locate and return a given colormap in a colortable using a min-match
 * search against a given name.
 *
 * Input:
 *  ctab   Ctable *   The color table to look up the color map in.
 *  name     char *   The name to look up the color table with.
 * Output:
 *  return   Cmap *   The located colormap, or NULL if not found.
 */
Cmap *get_Cmap(Ctable *ctab, char *name)
{
  Cmap *cmap = (Cmap *) get_symbol(ctab->symtab, name, 1);
/*
 * If the colormap was found and is a colormap of class CM_COLOR,
 * make this colormap to default colormap.
 */
  if(cmap) {
    if(cmap->class == CM_COLOR)
      add_Cmap(ctab, "color", cmap);
    ctab->cmap = cmap;
  };
  return cmap;
}

/*.......................................................................
 * Install a new color map for subsequent use by PGIMAG.
 *
 * Input:
 *  cmap     Cmap *  The color map to be installed.
 *  contra  float    Contrast of color ramp (normally 1.0).
 *  bright  float    Brightness at the center color index (normally 0.5).
 * Output:
 *  return    int    The number of colors defined.
 */
int recolor(Cmap *cmap, float contra, float bright)
{
/*
 * Install the color map.
 */
  return plcmap(cmap->l, cmap->r, cmap->g, cmap->b, cmap->nc, contra, bright);
}

/*.......................................................................
 * Install a color table to be used by PGIMAG.
 *
 * Input:
 *  l      float *  An array of 'nc' normalized ramp levels.
 *                  Levels must be sorted in increasing order.
 *                   0.0 places a color at the beginning of the ramp.
 *                   1.0 places a color at the end of the ramp.
 *                  Colors outside these limits will not appear if
 *                  contra=1.0 and bright=0.5.
 *  r      float *  An array of 'nc' normalized red intensities.
 *  g      float *  An array of 'nc' normalized green intensities.
 *  b      float *  An array of 'nc' normalized blue intensities.
 *  nc       int    The number of color table entries.
 *  contra float    Contrast of color ramp (normally 1.0).
 *  bright float    Brightness at the center color index (normally 0.5).
 * Output:
 *  return   int    The number of color indexes defined.
 */
static int plcmap(float *l, float *r, float *g, float *b, int nc,
		  float contra, float bright)
{
  const float minctr = 1.0/256; /* Minimum absolute contrast */
  int minind;  /* Start color index */
  int maxind;  /* End color index */
  int ci;      /* The color index being processed */
  int ntotal;  /* The number of color indexes to be processed. */
  int nspan;   /* The number of color indexes spanned by the color table. */
  int below;   /* The color table index of the level below the required level */
  int above;   /* The color table index of the level above the required level */
  int forward; /* True is the color table is to be traverse from 0 to nc-1 */
  float ca,cb; /* Fractional range of color indexes covered by color table */
  float stretch; /* fabs(cb-ca) */
/*
 * Check arguments.
 */
  if(nc <= 0 || r==NULL || g==NULL || b==NULL) {
    lprintf(stderr, "plcmap: Too few colormap levels.\n");
    return 0;
  };
/*
 * Determine the range of color indexes to be used.
 */
  cpgqcir(&minind, &maxind);
/*
 * Count the number of color indexes to be processed.
 */
  ntotal = maxind - minind + 1;
/*
 * No definable colors?
 */
  if(ntotal < 1 || minind < 1)
    return 0;
/*
 * Convert from contrast to the stretch of the ramp.
 * Limit the contrast to prevent a divide-by zero error.
 */
  stretch = 1.0f / (fabs(contra) < minctr ? minctr : fabs(contra));
/*
 * The brightness is only defined between 0 and 1.
 */
  if(bright < 0.0) bright = 0.0;
  if(bright > 1.0) bright = 1.0;
/*
 * Convert from brightness and contrast to normalized color index
 * coordinates ca and cb, at which to place the start and end of the
 * color table.
 */
  if(contra >= 0) {
    ca = 1.0 - bright * (1.0 + stretch);
    cb = ca + stretch;
  } else {
    ca = bright * (1.0 + stretch);
    cb = ca - stretch;
  };
/*
 * Determine the number of color indexes spanned by the color table.
 */
  nspan = floor(fabs(cb-ca) * ntotal);
/*
 * Determine the direction in which the color table is to be traversed.
 */
  forward = ca <= cb;
/*
 * Initialize the indexes at which to start searching the color table.
 */
  below = nc-1;   /* Start index for traversing the table from nc-1 to 0. */
  above = 0;      /* Start index for traversing the table from 0 to nc-1. */
/*
 * Buffer PGPLOT commands until the color map has been completely installed.
 */
  cpgbbuf();
/*
 * Linearly interpolate color table RGB values onto each color index.
 */
  for(ci=minind; ci<=maxind; ci++) {
/*
 * Turn the color index into a fraction of the range minind -> maxind.
 */
    float ci_frac = (float)(ci-minind) / (float)(maxind-minind);
/*
 * Determine the color table position that corresponds to color index, ci.
 */
    float level = nspan>0 ? (ci_frac-ca) / (cb-ca) : (ci_frac<=ca ? 0.0 : 1.0);
/*
 * Determine the indexes of the two color table entries that straddle
 * this position.
 */
    if(forward) {
      while(above<nc && l[above]<level)
	above++;
      below = above-1;
    } else {
      while(below>=0 && l[below]>level)
	below--;
      above = below + 1;
    };
/*
 * If the indexes lie outside the table, use the color of the
 * nearest edge of the table.
 */
    if(below < 0) {
      level = 0.0f;
      below = above = 0;
    } else if(above > nc-1) {
      level = 1.0f;
      below = above = nc-1;
    };
/*
 * Linearly interpolate the primary color intensities from color table
 * entries 'below' and 'above'.
 */
    {
      float lwid = l[above] - l[below];
      float lpos = lwid > minctr ? (level - l[below]) / lwid : 0.0f;
      float red   = r[below] + (r[above] - r[below]) * lpos;
      float green = g[below] + (g[above] - g[below]) * lpos;
      float blue  = b[below] + (b[above] - b[below]) * lpos;
/*
 * Intensities are only defined between 0 and 1.
 */
      if(red < 0.0f) red = 0.0f;
      if(red > 1.0f) red = 1.0f;
      if(green < 0.0f) green = 0.0f;
      if(green > 1.0f) green = 1.0f;
      if(blue < 0.0f) blue = 0.0f;
      if(blue > 1.0f) blue = 1.0f;
/*
 * Install the new color representation.
 */
      cpgscr(ci, red, green, blue);
    };
  };
/*
 * Reveal the changed color map.
 */
  cpgebuf();
/*
 * Return the total number of color indexes processed.
 */
  return ntotal;
}

/*.......................................................................
 * Lookup a given transfer function type.
 *
 * Input:
 *  name     char *  The name of a transfer function.
 * Output:
 *  return Cmtran    The colormap transfer function type that corresponds
 *                   to 'name'.
 */
Cmtran get_Cmtran(char *name)
{
  static Symtab *tab=NULL;  /* The symbol table of transfer functions */
  int i;
/*
 * If the transfer function descriptors have not yet been installed in
 * a symbol table, install them now.
 */
  if(tab==NULL) {
    tab = new_Symtab(n_std_trans, "Transfer function", amb_report,
		     (SYM_DEL(*)) 0);
    for(i=0; tab && i<n_std_trans; i++) {
      if(add_symbol(tab, std_trans[i].name, (void *) &std_trans[i].type, 0) && i==0)
	tab = del_Symtab(tab);
    };
  };
  if(tab) {
    Cmtran *tran = (Cmtran *) get_symbol(tab, name, 1);
    if(tran)
      return *tran;
  };
/*
 * If not found, return a linear transfer function.
 */
  return TR_LINEAR;
}

/*.......................................................................
 * Look up the name that goes with a given transfer function enumerator.
 *
 * Input:
 *  tran    Cmtran   The enumerator to look up.
 * Output:
 *  return    char * The name of the transfer function referenced by
 *                   tran.
 */
char *name_Cmtran(Cmtran tran)
{
  int i;
  for(i=0; i<n_std_trans; i++) {
    if(std_trans[i].type == tran)
      return std_trans[i].name;
  };
  return "none";
}
