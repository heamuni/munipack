/*

  fitshead - list, add, edit and remove of records in FITS header

*/

#include <fitsio.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <string.h>

void help()
{
  fprintf(stdout,"fitshead - help\n");
}

void version()
{
  fprintf(stdout,"version\n");
}

int trecognize(char *value)
{
  double x;
  char *l;

  if( strcmp(value,"T") == 0 || strcmp(value,"F") == 0 )
    return TLOGICAL;

  if( sscanf(value,"%lf",&x) == 1 ) {

    for(l = value; *l != '\0'; l++ )
      if( *l == '.' || *l == 'e' || *l == 'E' )
	return TDOUBLE;
    
    return TINT;
  }

  return TSTRING;
}


int main(int argc, char *argv[])
{
  int opt,i,ninc;
  int op = 0;
  size_t nlen;
  char *keys = NULL;
  char *k,*l,**inclist;
  char *fname = NULL;
  char *value = NULL;
  char *comment = NULL;
  fitsfile *f;
  char record[FLEN_CARD];
  int status = 0, nhead;
  int dtype;
  int bval,ival;
  double dval;
  void *val;

  while ((opt = getopt(argc, argv, "k:a:e:l:dc:w:hV")) != -1) {
    fprintf(stderr,"%c \n",opt);
    switch (opt) {
    case 'a':
    case 'e':
    case 'w':
    case 'l':
      op = opt;
      value = optarg;
      //      flags = 1;
      break;
    case 'd':
      op = opt;
      break;      
    case 'c':
      comment = optarg;
      //      flags = 1;
      break;
    case 'k':
      keys = optarg;
      break;
    case 'h':
      help();
      return 0;
    case 'V':
      version();
      return 0;
    default: /* '?' */
      fname = optarg;
      /*
      fprintf(stderr, "Usage: %s [-t nsecs] [-n] name\n",
	      argv[0]);
      exit(EXIT_FAILURE);
      */
    }
  }

  for(i = optind; i < argc; i++)
    fname = argv[i];

  if( ! fname ) {
    fprintf(stderr,"Please, specify a FITS filename.\n");
    return 1;
  }

  //  fprintf(stderr,"%s %s\n",keys,fname);


  if( op == 0 ) {
    /* list mode */

    if ( !fits_open_file(&f, fname, READONLY, &status) ) {

      if( keys ) {
	/* print selected keys */
	
	/* compute number of ',' */
	ninc = 0;
	for(k = keys; *k != '\0'; k++)
	  ninc += *k == ',' ? 1 : 0;
	ninc++;

	inclist = malloc(sizeof(char *)*ninc);

	k = keys;
	for(i = 0; i < ninc - 1; i++) {
	  l = strchr(k,',');
	  nlen = l - k;
	  inclist[i] = malloc(sizeof(char)*(nlen+1));
	  strncpy(inclist[i],k,nlen);
	  k = l + 1;
	}
	inclist[ninc-1] = strdup(k);

	fits_read_record(f,0,record,&status);
	while( fits_find_nextkey(f,inclist,ninc,NULL,0,record,&status) 
	       != KEY_NO_EXIST && status == 0 ){
	  fprintf(stdout,"%s\n",record);
	}
	status = 0;

	for(i = 0; i < ninc; i++) 
	  free(inclist[i]);
	free(inclist);

      }
      else {
	/* print full header */

	fits_get_hdrspace(f,&nhead,NULL,&status);
	for(i = 0; status == 0 && i < nhead; i++) {
	  if( fits_read_record(f,i+1,record,&status) == 0 )
	    fprintf(stdout,"%s\n",record);
	}
      }

      fits_close_file(f, &status);
    }

    if (status) fits_report_error(stderr, status);
    return(status);

  }
  else {
    /* modify mode */

    if ( !fits_open_file(&f, fname, READWRITE , &status) ) {

      fprintf(stdout,"%s\n",value);

      if( keys ) {
	
	if( op == 'a' && strstr(keys,"COMMENT") )
	  fits_write_comment(f,value,&status);

	if( op == 'a' && strstr(keys,"HISTORY") )
	  fits_write_history(f,value,&status);

	else {

	  if ( value ) {

	    dtype = trecognize(value);
	    if( dtype == TLOGICAL ) {
	      bval = strcmp(value,"T") == 0 ? TRUE : FALSE;
	      val = &bval;
	    }
	    else if( dtype == TINT ) {
	      sscanf(value,"%d",&ival);
	      val = (void *) &ival;
	    }
	    else if( dtype == TDOUBLE ) {
	      sscanf(value,"%lf",&dval);
	      val = &dval;
	    }
	    else {
	      dtype = TSTRING;
	      val = value;
	    }
	  }
	  printf("%d\n",dtype);

	  //	  trecognize(value,&dtype,&nhead);
	  //	  val = &ival;
	  //	printf("%d '%s' %d %d %d\n",dtype,value,*(int *)val,nhead,ival);

	  if( op == 'a' )
	    fits_write_key(f,dtype,keys,val,comment,&status);
	  if( op == 'e' )
	    fits_update_key(f,dtype,keys,val,comment,&status);
	  if( op == 'd' )
	    fits_delete_key(f,keys,&status);
	  if( op == 'l' )
	    fits_modify_name(f,keys,value,&status);
	}

      }
      else {

	if( op == 'w' ) {
	  if( strcmp(value,"-") == 0 ) {
	    while( fgets(&record[0],FLEN_CARD,stdin) ) {
	      for(i = 0; i < strlen(record); i++)
	      /*	      for(i = 0; i < FLEN_CARD; i++)*/
		/*		if( iscntrl(record[i]) )*/
		if( record[i] == '\n' )
		  record[i] = ' ';
	      record[FLEN_CARD-1] = '\0';
	      fits_write_record(f,record,&status);
	    }
	  }
	  else
	    fits_write_record(f,value,&status);
	}
      }

      fits_close_file(f, &status);
    }

    if (status) fits_report_error(stderr, status);
    return(status);

  }
  return 0;
}

