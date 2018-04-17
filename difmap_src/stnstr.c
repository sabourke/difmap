#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "logio.h"
#include "obs.h"
#include "vlbutil.h"

typedef struct {
  char *station;   /* Prefix name to recognise station by */
  char *array;     /* If the station prefix is redundant name the array */
  char *abrev;     /* The abbreviation to refer to the telescope by */
} Sabrev;

static Sabrev *stnabr(char *name);

static Sabrev unknown = {"Unknown", NULL, "?"};
 
static Sabrev stntab[]={
  {"HRAS",  NULL,   "F"},{"GRAS",  NULL,   "F"},{"NRAO",  NULL,   "G"},
  {"HAY",   NULL,   "K"},{"HST",   NULL,   "K"},{"ILL",   NULL,   "V"},
  {"EFF",   NULL,   "B"},{"ALG",   NULL,   "C"},{"ARO",   NULL,   "C"},
  {"PU",    NULL,  "Pu"},{"CRI",   NULL,   "R"},{"KRI",   NULL,   "R"},
  {"SIM",   NULL,   "R"},{"ONS",   NULL,   "S"},{"CHI",   NULL,   "U"},
  {"DWI",   NULL,   "W"},{"SAF",   NULL,   "E"},{"HAR",   NULL,   "E"},
  {"WESTF", NULL,   "K"},{"MARY",  NULL,   "N"},{"VLA",   "VLA",  "Y"},
  {"WSRT0", NULL,  "W0"},{"MED",   NULL,   "L"},{"BGNA",  NULL,   "L"},
  {"BOL",   NULL,   "L"},{"NOB",   NULL,   "M"},{"NRO",   NULL,   "M"},
  {"JAP",   NULL,   "M"},{"NOTO",  NULL,  "No"},{"TOR",   NULL,   "Z"},
  {"DEF",   NULL,  "Df"},{"CAM",   NULL,  "Cb"},{"LOV",   NULL,  "J1"},
  {"JBNK1", NULL,  "J1"},{"JBNK2", NULL,  "J2"},{"MET",   NULL,   "V"},
  {"FIN",   NULL,   "V"},{"ITA",   NULL,   "X"},{"ATI",   NULL,   "X"},
  {"OOTY",  NULL,  "Oo"},{"SHA",   NULL,  "Sh"},{"KASH",  NULL,  "Ka"},
  {"PER",   NULL,  "Pr"},{"ALI",   NULL,  "As"},{"PAR",   NULL,  "Pk"},
  {"PKS",   NULL,  "Pk"},{"CUL",   NULL,  "Cg"},{"HOB",   NULL,  "Hb"},
  {"NAN",   NULL,  "Nc"},{"MAD",   NULL,   "D"},{"DSS13", NULL,  "Dv"},
  {"DSS14", NULL,  "Dm"},{"DSS15", NULL,  "Dg"},{"DSS6",  NULL,  "Ds"},
  {"DSS4",  NULL,  "Dt"},{"PIE",   NULL,  "Pt"},{"KIT",   NULL,  "Kp"},
  {"LOS",   NULL,  "La"},{"VLBA_PT", NULL,"Pt"},{"VLBA_KP", NULL,"Kp"},
  {"VLBA_LA", NULL,"La"},{"VLBA_FD", NULL,"Fd"},{"VLBA_NL", NULL,"Nl"},
  {"VLBA_BR", NULL,"Br"},{"VLBA_OV", NULL,"Ov"},{"VLBA_SC", NULL,"Sc"},
  {"VLBA_HN", NULL,"Hn"},{"VLBA_MK", NULL,"Mk"},{"BONN", NULL,    "B"},
  {"OVRO",  NULL,   "O"},{"AN", "mma", "Mm"},{"WSRT",  "WSRT",   "W"},
  {0,0}};

static char *subarray_string(Subarray *sub, char *subnam, int slen);

/*.......................................................................
 * Given a station name return its standard abbreviation.
 *
 * Input:
 *   name     char *  The station name.
 * Output:
 *   return Sabrev *  The container of the official telescope name and
 *                    its abbreviation.
 *                    If the telescope is not recognised &unkown will
 *                    be returned.
 */
static Sabrev *stnabr(char *name)
{
  static char wrkstr[10];
  static char *inptr;
  static char *ouptr;
  static Sabrev *sptr;
  static int i;
/*
 * Make an upper case copy of the input name.
 */
  inptr=name;
  ouptr= &wrkstr[0];
  for(i=0; *inptr != '\0' && i<sizeof(wrkstr)-1; i++)
    *(ouptr++) = (char) toupper((unsigned int) *(inptr++));
  *ouptr = '\0';
/*
 * Locate the abbrevation.
 */
  for(sptr = &stntab[0]; sptr->station != 0; sptr++) {
    inptr=wrkstr;
    ouptr=sptr->station;
    for(;*inptr != '\0' && *ouptr != '\0';inptr++,ouptr++)
      if(*inptr != *ouptr)
	break;
    if(*ouptr == '\0')
      return sptr;
  };
/*
 * Station not found.
 */
  return &unknown;
}

/*.......................................................................
 * Create a string naming each subarray.
 *
 * Input:
 *  ob  Observation *  The observation whose stations are to be listed.
 *  arrstr     char *  The return string of length 'slen'.
 *  slen        int    The available length of string 'arrstr'.
 * Output:
 *  return      int    0 - normal completion.
 *                     1 - string too short for station list.
 */
int stnstr(Observation *ob, char *arrstr, int slen)
{
  char *last=NULL; /* The pointer to the last subarray name in 'arrstr' */
  char next[40];   /* The name of the next subarray */
  char *aptr;      /* A pointer into 'arrstr' */
  int nleft=slen-1;/* Number of unused chars in 'arrstr' not including '\0' */
  int i;
/*
 * Check arguments.
 */
  if(!ob_ready(ob, OB_INDEX, "stnstr"))
    return 1;
/*
 * Get the abbreviation for each telescope in each sub-array of 'ob'.
 */
  aptr = arrstr;
  for(i=0; i<ob->nsub; i++) {
    Subarray *sub = ob->sub + i;
    int sublen;
/*
 * Get the next array name.
 */
    if(subarray_string(sub, next, sizeof(next)) == NULL)
      return 1;
/*
 * If it is the same as the last one that was written into arrstr[],
 * don't add it redundantly.
 */
    if(last && strcmp(last, next) == 0)
      continue;
/*
 * Is there room for a space, the next sub-array name and a terminating
 * '\0' character?
 */
    sublen = strlen(next);
    if(sublen + 2 > nleft) {
      lprintf(stderr, "The array title has been truncated.\n");
      return 1;
    };
/*
 * Add a space between the sub-array names.
 */
    *aptr++ = ' ';
    nleft -= 1;
/*
 * Append the new subarray name.
 */
    strcpy(aptr, next);
/*
 * Record a pointer to the the newly copied substring for comparison
 * with the next array name, then prepare for the next subarray.
 */
    last = aptr;
    aptr += sublen;
    nleft -= sublen;
  };
  return 0;
}

/*.......................................................................
 * Compose a string containing the array name or VLBI-style telescope
 * list of a given sub-array.
 *
 * Input:
 *  sub    Subarray *  The subarray to be named.
 *  subnam     char *  The return string of length 'slen'.
 *  slen        int    The available length of string 'subnam'.
 * Output:
 *  return     char *  This will be 'subnam' unless there was an error,
 *                     in which case NULL will be returned.
 */
static char *subarray_string(Subarray *sub, char *subnam, int slen)
{
  char *arrnam;   /* The optional array name from a binary antenna table */
  int i;
/*
 * Was an array name specified in the antenna table?
 */
  arrnam = sub->binan ? sub->binan->arrnam : NULL;
/*
 * If the antenna table specified an array name, use it unless it is
 * VLBI or VLBA. In the latter cases the stations will be named individually,
 * as is the custom for these arrays.
 */
  if(arrnam && *arrnam &&
     strcmp(arrnam, "VLBI") != 0 && strcmp(arrnam, "VLBA") != 0) {
/*
 * Check wether the array name will fit into the output string.
 */
    if(strlen(arrnam) + 1 > slen) {
      lprintf(stderr, "stnstr: Displaying truncated station string.\n");
      return NULL;
    };
/*
 * Copy the array name into the output string.
 */
    strcpy(subnam, arrnam);
/*
 * If an array name was not given, or the array was VLBI or VLBA, we will
 * compose a string containing a concatenated list of the standard
 * abbreviations for the stations of the subarray.
 */
  } else {
    int nused = 0;      /* The number of used characters in subnam[] */
    int nleft = slen;   /* The number of characters free in subnam[] */
    Sabrev *first=NULL; /* The table entry in stntab[] of the first antenna */
    int onetel = 1;     /* True if all antennas have the same abbreviation */
/*
 * Append abbreviations for each telescope of the subarray, into subnam[].
 */
    for(i=0; i<sub->nstat && nleft>1; i++) {
/*
 * Get the abbreviation and its length.
 */
      Sabrev *sptr = stnabr(sub->tel[i].name);
      int alen = strlen(sptr->abrev);
/*
 * If no abbreviation is known, we will substitute the first letter
 * of the telescope name.
 */
      if(sptr == &unknown)
	alen = 1;
/*
 * Will it fit in the output string?
 */
      if(alen + 1 > nleft) {
	lprintf(stderr, "The array title has been truncated.\n");
	return NULL;
      };
/*
 * Append the abbreviation to the output string, using the first character
 * of the telescope name if no abbreviation was found.
 */
      if(sptr==&unknown) {
	subnam[nused] = toupper((unsigned int) sub->tel[i].name[0]);
	subnam[nused+1] = '\0';
      } else {
	strcpy(subnam + nused, sptr->abrev);
      };
      nused += alen;
      nleft -= alen;
/*
 * Keep a record of whether any of the abbreviations differ. If they
 * don't then the abbreviation is probably a composite array name.
 */
      if(first==NULL)
	first = sptr;
      onetel = onetel && sptr != &unknown && sptr == first;
    };
/*
 * Did all the stations have the same recognised abbreviation, and is this
 * identified as an array?
 */
    if(onetel && first->array != NULL) {
      int alen = strlen(first->array);
/*
 * Will the station name fit in the output string?
 */
      if(alen + 1 > slen) {
	lprintf(stderr, "The array title has been truncated.\n");
	return NULL;
      };
/*
 * Overwrite the station abbreviations with the official array name.
 */
      strcpy(subnam, first->array);
    };
  };
  return subnam;
}
