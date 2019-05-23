
// file.dicz.R contains rules for dictionary words, first int is # of terminals
// file.dicz.C contains the sequences, separated by terminals >= 256
// for each sequence, create balanced rules to convert it into one nonterminal
// the terminals in file.parse must be renamed as these one nonterminals
// and the nonterminals must be shifted accordingly
// the rules of .parse.R and .dicz.R must be appended

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>

int bits (int x)

   { int l=0;
     while (x) { x>>=1; l++; }
     return l;
   }

void main (int argc, char **argv)

  { FILE *diczR,*diczC,*parseR,*parseC,*R,*C;
    char fname[1024];
    int terms; // how many terminals in diczR
    int rules; // how many rules in diczR
    int prules; // how many rules in parseR
    int sizeC; // size of diczC
    int psizeC; // size of parseC
    int phrases; // how many phrases in dicz = terminals in parseR
    struct stat s;
    int v256 = 256;
    int i,p;
    int val[2];
    int *AdiczC; // to read diczC in memory
    int *transl; // translation table for phrases of dicz to their nonterms

    puts("==== Command line:");
    for(i=0;i<argc;i++)
      printf(" %s",argv[i]);
    puts("\n");

    if (argc != 2)
       { fprintf (stderr,
		  "Usage: %s <file> makes <file>.[RC] from <file>.[dicz+parse].[RC]\n",
		  argv[0]);
 	 exit(1);
       }

	// open all *R files and read basic data from them

    sprintf(fname,"%s.R",argv[1]);
    R = fopen(fname,"w");
    if (R == NULL)
       { fprintf (stderr, "Cannot open %s\n",fname);
 	 exit(1);
       }

    sprintf(fname,"%s.dicz.int.R",argv[1]);
    diczR = fopen(fname,"r");
    if (diczR == NULL)
       { fprintf (stderr, "Cannot open %s\n",fname);
 	 exit(1);
       }
    fread(&terms,sizeof(int),1,diczR);
    stat (fname,&s);
    rules = (s.st_size-sizeof(int))/(2*sizeof(int));

    sprintf(fname,"%s.parse.R",argv[1]);
    parseR = fopen(fname,"r");
    if (parseR == NULL)
       { fprintf (stderr, "Cannot open %s\n",fname);
 	 exit(1);
       }
    fread (&phrases,sizeof(int),1,parseR); // must hold phrases + 256 = terms
    stat (fname,&s);
    prules = (s.st_size-sizeof(int))/(2*sizeof(int));

	// shift all rules of diczR, removing terminals in 256..terms-1
	// leaves cursor at the end of R, ready to add more rules

    transl = (int*) malloc ((phrases+1)*sizeof(int));
    if ((transl == NULL) && (phrases > 0))
       { fprintf (stderr, "Cannot allocate %li bytes for the phrases of %s.dicz.int.R\n",
		  phrases*sizeof(int),argv[1]);
 	 exit(1);
       }

    fwrite(&v256,sizeof(int),1,R); // could be better, we assume 256 terminals 0..255
    for (i=0;i<rules;i++)
	{ fread(val,sizeof(int),2,diczR);
 	  if (val[0] >= 256) val[0] -= terms-256; 
 	  if (val[1] >= 256) val[1] -= terms-256; 
	  fwrite(val,sizeof(int),2,R);
	}
    fclose(diczR);

	// read diczC in memory, identify the terms and create balanced rules on them
	// adding those to R and remembering the mapping phrase -> top nonterm in transl

    sprintf(fname,"%s.dicz.int.C",argv[1]);
    diczC = fopen(fname,"r");
    if (diczC == NULL)
       { fprintf (stderr, "Cannot open %s\n",fname);
 	 exit(1);
       }
    stat(fname,&s);
    sizeC = s.st_size/sizeof(int);
    AdiczC = (int*) malloc (s.st_size);
    if (AdiczC == NULL)
       { fprintf (stderr, "Cannot allocate %li bytes for %s.dicz.int.C\n",s.st_size,argv[1]);
 	 exit(1);
       }
    fread (AdiczC,sizeof(int),sizeC,diczC);
    fclose (diczC);

    i = 0; p = 1;
    while (i < sizeC)
       { int j = i;
	 int ni;
	 while ((AdiczC[j] < 256) || (AdiczC[j] >= terms))
	    { if (AdiczC[j] >= 256) AdiczC[j] -= terms-256; 
	      j++;
	    }
	 ni = j+1;
	 while (j-i > 1)
	   { int k = i;
	     int ko = i;
	     while (k+1 < j)
		{ val[0] = AdiczC[k++]; val[1] = AdiczC[k++];
		  fwrite(val,sizeof(int),2,R);
		  AdiczC[ko++] = 256 + rules++;
		}
	     if (k < j) AdiczC[ko++] = AdiczC[k];
	     j = ko;
	   } 
	 transl[p++] = AdiczC[i];
	 i = ni;
       }
    free (AdiczC);

	// translate the rules in .parse.R according to table transl

    for (i=0;i<prules;i++)
       { fread(val,sizeof(int),2,parseR);
	 if (val[0] < phrases) val[0] = transl[val[0]];
	 else val[0] = val[0] - phrases + rules + 256;
	 if (val[1] < phrases) val[1] = transl[val[1]];
	 else val[1] = val[1] - phrases + rules + 256;
	 fwrite(val,sizeof(int),2,R);
       }
    fclose (parseR);
    if (fclose(R) != 0)
       { fprintf (stderr,"Cannot close %s.R\n",argv[1]);
	 exit(1);
       }
    
	// translate the rules in .parse.C according to table transl, open C files

    sprintf(fname,"%s.parse.C",argv[1]);
    parseC = fopen(fname,"r");
    if (parseC == NULL)
       { fprintf (stderr, "Cannot open %s\n",fname);
 	 exit(1);
       }
    stat (fname,&s);
    psizeC = s.st_size/sizeof(int);

    sprintf(fname,"%s.C",argv[1]);
    C = fopen(fname,"w");
    if (C == NULL)
       { fprintf (stderr, "Cannot open %s\n",fname);
 	 exit(1);
       }

    for (i=0;i<psizeC;i++)
       { fread(val,sizeof(int),1,parseC);
	 if (val[0] < phrases) val[0] = transl[val[0]];
	 else val[0] = val[0] - phrases + rules + 256;
	 fwrite(val,sizeof(int),1,C);
       }
    free (transl);
    fclose (parseC);
    if (fclose(C) != 0)
       { fprintf (stderr,"Cannot close %s.C\n",argv[1]);
	 exit(1);
       }

    stat(argv[1],&s);
    rules += prules;
    fprintf(stderr,"Prefix-Free + Repair succeeded\n\n");
    fprintf(stderr,"  Original file size: %li\n",s.st_size);
    fprintf(stderr,"  Number of rules: %i\n",rules);
    fprintf(stderr,"  Final sequence length: %i\n",psizeC);
    fprintf(stderr,"  Estimated output size (bytes): %ld\n",
      (long) ( (4.0*rules+((float)bits(256+rules))*(rules+psizeC))/8 )+1);
    fprintf(stderr,"  Compression ratio: %0.2f\%\n",
                    100.0*(4.0*rules+((float)bits(256+rules))*(rules+psizeC))
                        /(8.0*s.st_size));

    puts("=== postprocessing completed!");

    exit(0);
  }
