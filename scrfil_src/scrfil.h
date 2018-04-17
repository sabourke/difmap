#ifndef scrfil_h
#define scrfil_h

/* Scratch file utilities */

/* Return version numbered file name */

char *scrname(const char *base);

/* True if named file exists and is readable */

int file_exists(const char *name);

/* Invoke an external editor on a given file */

int ed_file(const char *name);

#endif
