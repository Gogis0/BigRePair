
#include <stdlib.h>
#include <stdio.h>

int main (int argc, char **argv)

  { char foname[1024];
    FILE *fi,*fo;
    int n=256;
    int c,i;

    puts("==== Command line:");
    for(i=0;i<argc;i++)
      printf(" %s",argv[i]);
    puts("\n");

    fi = fopen(argv[1],"r");
    if (fi == NULL)
  { fprintf(stderr,"Cannot open file %s\n",argv[1]);
    exit(1);
  }
    sprintf(foname,"%s.int",argv[1]);
    fo = fopen(foname,"w");
    if (fo == NULL)
  { fprintf(stderr,"Cannot create file %s\n",foname);
    exit(1);
  }
    while (1)
      { c = getc(fi);
  if ((c==0) || (c==-1)) break;
    if (c == 1) c = n++;
  fwrite (&c,sizeof(int),1,fo);
      }
    fclose(fo);
    fclose(fi);
    printf ("%i strings\n",n-256);
    puts("=== Preprocessing completed!");
    exit(0);
  }
