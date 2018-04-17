#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "utils.h"
/*
 * Local function protypes:
 */
static int match_single(char *start_ptr, char *end_ptr, char ch,
			 int *was_error);
static char get_next_char(char *string, char **tail, int *was_escaped);

/*.......................................................................
 * This is the top level matching routine. It returns 1 if the string
 * regexp matches the string, string - otherwise it returns 0.
 * The final argument is an error flag. It is normally 0 but when there
 * has been a syntax error it becomes 1.
 */
int match(char *regexp, char *string, int *was_error)
{
  int was_escaped;  /* Returned by get_next_char() to report that the */
                    /* character that it returned was escaped */
  char *start_ptr;  /* Start pointer of regexp for match_single() */
  char *end_ptr;    /* End pointer of regexp for match_single() */
  char *tail;       /* Pointer to tail part of some string not yet processed */
  char new_char;    /* Current char returned by get_next_char() */
  int match_test;   /* 1 or 0 depending on whether success is measured */
                    /* by a match or by no match */
/*
 * No error yet.
 */
  *was_error=0;
/*
 * Parse regexp.
 */
  for(;;) {
    new_char=get_next_char(regexp, &regexp, &was_escaped);
/*
 * Interpret new regexp character.
 */
    if(!was_escaped) {
      switch (new_char) {
/*
 * Stop at end of regexp string. If the string is also exhausted then
 * the match has succeded.
 */
      case '\0':
	return (*string == '\0') ? 1:0;
	break;
/*
 * Start of a single character regexp.
 */
      case '[':
/*
 * Keep a record of where the inside of the [] regexp started.
 */
	start_ptr = regexp;
/*
 * Locate the matching (un-escaped) ']'.
 */
	for(;;) {
	  new_char=get_next_char(regexp, &regexp, &was_escaped);
	  if(new_char== '\0' || (new_char==']' && !was_escaped) )
	    break;
	};
/*
 * If not found, report error and return.
 */
	if(new_char == '\0') {
	  lprintf(stderr, "Syntax error: \'[\' not matched in regexp\n");
	  *was_error=1;
	  return 0;
	};
/*
 * Record the pointer to the last char in the [] regexp.
 */
	end_ptr = regexp-2;
/*
 * If the first character of the regexp is an un-escaped '^' then
 * the result of the one-char regexp match is to be complemented.
 */
	new_char=get_next_char(start_ptr, &tail, &was_escaped);
	if(new_char == '^' && !was_escaped) {
	  start_ptr=tail;
	  match_test=1;
	}
	else {
	  match_test=0;
	};
/*
 * Make sure that the [] expression contains something.
 */
	if(end_ptr < start_ptr) {
	  lprintf(stderr, "Syntax error: Empty [] regexp encounterred\n");
	  *was_error=1;
	  return 0;
	};
/*
 * If the character following the closing ']' is an unescaped '*', then
 * the single-character regexp must match 0 or more characters. Otherwise
 * it should match exactly one.
 */
	new_char=get_next_char(regexp, &tail, &was_escaped);
	if(new_char=='*' && !was_escaped) {
	  regexp=tail;
	  for(;;) {
	    if(match(regexp, string, was_error) == 1)
	      return 1;
	    if(*was_error) return 0;
	    if(*string == '\0')
	      return 0;
	    if(match_single(start_ptr, end_ptr, *(string++), was_error)
	       == match_test)
	      return 0;
	  };
	}
	else {
/*
 * A '+' following a [] expression says that it should match exactly
 * one character. This is only really required when one wants to
 * a following '*' not to be associated with the [].
 */
	  if(new_char=='+' && !was_escaped) regexp=tail;
	  if(*string == '\0')
	    return 0;
	  if(match_single(start_ptr, end_ptr, *(string++), was_error)
	     == match_test)
	    return 0;
	};
	continue;
/*
 * '.' matches any character in string.
 */
      case '.':
	if(*string == '\0')
	  return 0;
        string++;
	continue;
/*
 * Zero or more characters in string.
 */
      case '*':
/*
 * Further '*'s are redundant.
 */
	while(get_next_char(regexp, &tail, &was_escaped)=='*'
	      && !was_escaped) regexp=tail;
/*
 * End of regexp string? If so the rest of string definitely matches.
 */
	if(*regexp == '\0') return 1;
/*
 * Now we don't know where the next part of the regexp will continue
 * in string, so keep recursively calling match() for each subsequent
 * character in string until we hit the end of string (in which case
 * the match has failed) or match() returns success.
 */
	while(*string != '\0') {
	  if(match(regexp, string++, was_error)==1)
	    return 1;
	};
	return 0;
/*
 * The character is not a special one - nor are the escaped characters
 * that didn't get into this switch, they will all be handled together
 * off the bottom of the switch expression.
 */
      default:
	break;
      };
    };
/*
 * Simple character in regular expression must match that in string.
 * If the character is a \ however it is assumed to escape the following
 * character - which in turn must match the character in the string.
 */
    if(*string == new_char) {
      string++;
    }
    else {
      return 0;
    };
  };
}

/*.......................................................................
 * Given an input string, return either the first character or if that
 * first character is a \ return the corresponding escaped character.
 * If there is an error return '\0'. The pointer into the input string
 * will be returned via the second argument, incremented to the point
 * just after the last character used. The flag, 'was_escaped' is returned
 * as 1 if the returned character was an escape sequence, and 0 otherwise.
 */
static char get_next_char(char *string, char **tail, int *was_escaped)
{
  char ch;   /* Will hold the character that is to be returned */
/*
 * Check for escape sequences.
 */
  switch (*string) {
  case '\\':
    string++;
    *was_escaped=1;
    switch (*string) {
/*
 * Standard escape equivalents for control characters.
 */
    case 'n':
      ch = '\n';
      break;
    case 'r':
      ch = '\r';
      break;
    case 't':
      ch = '\t';
      break;
    case 'f':
      ch = '\f';
      break;
/*
 * Octal escape sequence.
 */
    case '0': case '1': case '2': case '3': case '4': case '5': case '6':
    case '7': case '8': case '9':
      ch = (char) strtol(string, &string, 8);
      break;
/*
 * Hex escape sequence.
 */
    case 'x': case 'X':
      if(!isxdigit((int) string[1])) {
	lprintf(stderr, "Incomplete \\x.. hexadecimal escape sequence\n");
	*tail=string;
	return '\0';
      };
      ch = (char) strtol(string, &string, 16);
      break;
/*
 * No char to escape - whoops - error.
 */
    case '\0':
      lprintf(stderr, "Incomplete escape sequence at end of string\n");
      *tail = string;
      return '\0';
      break;
/*
 * Any other character gets passed unscathed.
 */
    default:
      ch = *string;
    };
    break;
  default:
    ch = *string;
    *was_escaped=0;
    break;
  };
/*
 * Return the result.
 */
  *tail = ++string;
  return ch;
}

/*.......................................................................
 * Try to match a single character ch with a [] regexp character list.
 * Return 1 on success, 0 on failure. Normally *was_error is returned
 * as 0, but if there is a regexp syntax error it is returned as 1.
 * The regexp can be formed of any arrangement of the following:
 *
 * A-Z  : Any capital letter between the two characters provided.
 * a-z  : Any lower case letter between the two characters provided.
 * 1-9  : Any digit between the numbers specified.
 * adbc : Any charcter from the specified list.
 * ^    : When placed before any of the above, the match will succede
 *      : if the character is not one of those specified. Characters
 *      : may be escaped to remove any special meanings or to include
 *      : control characters etc..
 */
static int match_single(char *start_ptr, char *end_ptr, char ch,
			int *was_error)
{
  int in_range;  /* Flags when half way through processing a regexp range */
  int was_escaped;/* Flags characters that were escaped when read */
  char last_char;/* Keeps a record of first char of a potential regexp range */
  char new_char; /* Latest char read from regexp */
  char *regexp;  /* Used to step through regexp */
  char *ptr;     /* Pointer for general usage */
/*
 * Strings holding sensible collating sequence for ranges.
 * It would be much faster to assume ASCII but not very portable.
 */
  static char upper[]="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  static char lower[]="abcdefghijklmnopqrstuvwxyz";
  static char digit[]="0123456789";
/*
 * No errors yet.
 */
  *was_error=0;
/*
 * Copy start pointer to regexp-stepping pointer.
 */
  regexp=start_ptr;
/*
 * Process each char list in regexp in turn until a match is found
 * or the ending ']' is reached.
 */
  in_range=0;
  last_char='\0';
  while(regexp <= end_ptr) {
/*
 * Get the next character of the regexp, with due regard for
 * escape sequences.
 */
    new_char=get_next_char(regexp, &regexp, &was_escaped);
/*
 * Check for the range-symbol.
 */
    if(!was_escaped && new_char == '-') {
      if(in_range || last_char == '\0') {
	lprintf(stderr, "Incomplete regexp range\n");
	*was_error=1;
	return 0;
      };
      in_range=1;
      continue;
    };
/*
 * Try to match the latest character with ch.
 */
    if(ch==new_char)
      return 1;
/*
 * Is this the final char of a range designation?
 */
    if(!in_range) {
      last_char=new_char;
    }
/*
 * Complete range designation recieved - process it.
 */
    else {
/*
 * Check types of the two characters and see if ch is between them.
 */
      if(isdigit((int) last_char) && isdigit((int) new_char)) {
	ptr = strchr(digit, ch);
	if(strchr(digit,last_char) <= ptr && strchr(digit, new_char) >= ptr)
	  return 1;
      }
      else if(isupper((int) last_char) && isupper((int) new_char)) {
	ptr = strchr(upper, ch);
	if(strchr(upper,last_char) <= ptr && strchr(upper, new_char) >= ptr)
	  return 1;
      }
      else if(islower((int) last_char) && islower((int) new_char)) {
	ptr = strchr(lower, ch);
	if(strchr(lower,last_char) <= ptr && strchr(lower, new_char) >= ptr)
	  return 1;
      }
      else {
	lprintf(stderr, "Syntax error in regexp character range\n");
	*was_error=1;
	return 0;
      };
/*
 * Range parsed, but with no match - prepare for next part of
 * regexp.
 */
      in_range = 0;
      last_char = '\0';
    };
  };
/*
 * No match found.
 */
  return 0;
}
