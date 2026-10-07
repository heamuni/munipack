
/*
  Compile as:

  $ g++ -Wall fortranio.cpp fi2.cpp

*/

#include "fortranio.h"
#include <iostream>
#include <string>
#include <vector>
#include <cstring>

using namespace std;



int main()
{
  char *a[12];
  //  a = malloc(10*sizeof(*char));
  a[0] = strdup("324.092 5.572    0 0.143");
  a[1] = strdup(" '''If kocka''s wants.' ");
  a[2] = strdup(" 'If kocka''s wants. ");
  a[3] = strdup(" 'If kocka''s wants.' ");
  a[4] = strdup("'''If kocka''s wants.'''''");
  a[5] = strdup("  ' If kocka'' wants. ' ");
  a[6] = strdup(" If kocka''s wants. ");
  a[7] = strdup("3 'If kocka''s wants.' 56.3,0.36 4,5 'aa'");
  a[8] = strdup("'TEST' 'FILTER' '''R''' ''");
  a[9] = strdup("'TRIMSEC' '[1:2,3:4]'");
  a[10] = strdup("''[1:2,3:4]''");
  a[11] = NULL;

  for(size_t l = 0; a[l] != NULL ; l++) {
    cout << a[l] << endl;
    vector<string> t(strftok(a[l]));
    for(size_t i = 0; i < t.size(); i++)
      cout << "#" << i << " " << t[i] << endl;
  }

  return 0;
}
