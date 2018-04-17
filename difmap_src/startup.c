#include <stdio.h>
#include <stdlib.h>
#include "logio.h"
#include "sphere.h"

/* List all modules to be installed in the symbol tables */

extern Module m_general;
extern Module m_graphics;
extern Module m_maths;
extern Module m_iolib;
extern Module m_difmap;

static Module *modules[] = {
  &m_general,
  &m_graphics,
  &m_maths,
  &m_iolib,
  &m_difmap
};

/*.......................................................................
 * Start difmap.
 */
int cdifmap_(void)
{
/*
 * Write start-up message.
 */
#include "version.h"
  lprintf(stdout, "Copyright (c) 1993-2008 California Institute of Technology. All Rights Reserved.\n");
  lprintf(stdout, "Type 'help difmap' to list difference mapping commands and help topics.\n");
/*
 * Start up the interface.
 */
  return startup(modules, sizeof(modules)/sizeof(Module *), "DIFMAP_LOGIN");
}

/*.......................................................................
 * This function gets called from f77main.f on systems that don't postfix
 * FORTRAN referenced functions with underscores. On most systems
 * cdifmap_() gets called directly from f77main.c and this dummy function
 * is not called.
 */
int cdifmap(void)
{
  return cdifmap_();
}
