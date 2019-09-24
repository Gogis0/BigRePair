#include <stdlib.h>
#include <stdio.h>

// transform a dictionary converting chars to int32_t's and adding a 
// a unique terminator after each string
// the values used for unique terminators are 256, 257, and so on ... 
// this version uses the old (Spire'19) dictionary format where 
// each dictionary string is terminated by 0x1 and there is a 0x0
// at the end of the dictionary   

int main (int argc, char **argv)
{ 
  char foname[1024];
  FILE *fi,*fo;
  int n=256;  // first integer value used as a separator
  int c,i;

  puts("==== Command line:");
  for(i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("\n");

  // open dictionary file
  fi = fopen(argv[1],"r");
  if (fi == NULL) { 
    fprintf(stderr,"Cannot open file %s\n",argv[1]);
    exit(1);
  }
  sprintf(foname,"%s.int",argv[1]);
  // open integer output file
  fo = fopen(foname,"w");
  if (fo == NULL) { 
    fprintf(stderr,"Cannot create file %s\n",foname);
    exit(1);
  }
  
  // main loop 
  while (1) { 
    c = getc(fi); // next char
    if ((c==0) || (c==-1)) break; // if 0 or EOF exit loop
    if (c == 1)   // if 1 then a dictionary string has ended
       c = n++;   // replace 1 with unique terminator 
    fwrite (&c,sizeof(int),1,fo); // write char or terminator as an int to output file 
  }
  fclose(fo);
  fclose(fi);
  printf ("%i strings\n",n-256);
  puts("=== Preprocessing completed!");
  exit(0);
}
