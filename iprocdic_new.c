#include <stdlib.h>
#include <stdio.h>


// iprocdic_new.c
//
// Intermediate tool for bigRepair for integers 
// derived from procdic_new.c and modified for the case
// in which the input alphabet of prefix free parsing consists of integers

// transform the PFP dictionary and adding a 
// a unique terminator after each string
// the values used for unique terminators are U, U+1, and so on ...
// where U is currently 1<<30 (ie 2^30).
// the program check that input symbols are indeed smaller than U 
// otherwise it exits with an error message
 
// This version uses the new (Spire'19) dictionary format where 
// the lengths (in symbols) of dictionary strings are provided in a 
// separate .len file in int32_t integers. 

#define Unique (1<<30)


int main (int argc, char **argv)
{ 
  char foname[1024], flname[1024];
  FILE *fi,*fo, *fl;
  int n=Unique;  // first integer value used as a separator
  int c,i;

  if(argc!=2) {
     printf("Usage: %s <dictionaryfilename>\n\n", argv[0]);
     exit(1);
  }

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
  // open dictionary lengths file
  sprintf(flname,"%s.len",argv[1]);
  fl = fopen(flname,"r");
  if (fl == NULL) { 
    fprintf(stderr,"Cannot open file %s\n",flname);
    exit(1);
  }
  // open integer output file
  sprintf(foname,"%s.int",argv[1]);
  fo = fopen(foname,"w");
  if (fo == NULL) { 
    fprintf(stderr,"Cannot create file %s\n",foname);
    exit(1);
  }
  
  // main loop 
  while (1) {
    // read length of next word 
    int wlen;
    int e = fread(&wlen,sizeof(int),1,fl);
    if(e!=1) {
      if(feof(fl)) break; // end file
      else {perror("Error reading dictionary word length"); exit(1); }
    }
    // read actual word
    for(int i=0;i<wlen;i++) {
      e = fread(&c,sizeof(int),1,fi);
      if(e!=1) {perror("Unexpected end of dictionary"); exit(1); }
      if(c<0 || c>=Unique) {fprintf(stderr, "Dictionary symbol %x larger than %x\n",c,Unique-1); exit(1);} 
      e = fwrite (&c,sizeof(int),1,fo); // write int to output file 
      if(e!=1) {perror("Error writing to .int file"); exit(1); }
    }
    e = fwrite (&n,sizeof(int),1,fo); // write terminator as an int to output file 
    if(e!=1) {perror("Error writing to .int file"); exit(1); }
    n++; // update terminator 
  }
  c = getc(fi); // there should be no further chars
  if(c!=EOF) {perror("Unexpected trailing chars in dictionary"); exit(1); }
  fclose(fo);
  fclose(fi);
  printf ("%i strings\n",n-Unique);
  puts("=== Preprocessing completed!");
  return 0;
}
