// ipostproc.c
//
// Intermediate tool for bigRepair for integers 
// derived from postproc.c and modified for the case
// in which the input alphabet of prefix free parsing consists of integers


// file.dicz.R contains rules for dictionary words, first int is # of terminals
// file.dicz.C contains the sequences, separated by terminals >= 256
// for each sequence, create balanced rules to convert it into one nonterminal
// the terminals in file.parse must be renamed as these one nonterminals
// and the nonterminals must be shifted accordingly
// the rules of .parse.R and .dicz.R must be appended

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <assert.h>

// constant larger than every true terminal
#define Unique (1<<30)

//  exit program with error msg if test is true
void die(int test,const char *msg)
{
  if(test) {
    perror(msg);
    exit(1);
  }
}


int bits (int x)

   { int l=0;
     while (x) { x>>=1; l++; }
     return l;
   }

// return the largest integer in f which is smaller than Unique
// error if a negative value is found 
int maxbelowUnique(FILE *f)
{
  int i,m=0;
  while(1) {
    int e = fread(&i,sizeof(int),1,f);
    if(e!=1) {
      die(ferror(f), "Read error (maxbelowUnique)");
      break;
    }  
    die(i<0,"Unexpected negative value (maxbelowUnique)");
    if(i<Unique && i>m)
      m=i;
  }
  return m;
}


int main (int argc, char **argv)
  { FILE *diczR,*diczC,*parseR,*parseC,*R,*C;
    char fname[1024];
    int terms; // how many terminals in diczR
    int rules; // how many rules in diczR
    int prules; // how many rules in parseR
    int sizeC; // size of diczC
    int psizeC; // size of parseC
    int phrases; // how many phrases in dicz = terminals in parseR
    struct stat s;
    int largest_symbol=0; // largest symbol in the original input alphabet
    int v256;   // legacy name for the size of the original input alphabet, ie largest_symbol+1  
    int i,p,e;
    int val[2];
    int *AdiczC; // to read diczC in memory
    int *transl; // translation table for phrases of dicz to their nonterms

    fprintf(stderr,"==== Command line:\n");
    for(i=0;i<argc;i++)
      fprintf(stderr," %s",argv[i]);
    fputs("\n",stderr);

    if (argc != 2)
       { fprintf (stderr,
      "Usage: %s <file>\nmakes <file>.[RC] from <file>.[dicz+parse].[RC]\n",
      argv[0]);
   exit(1);
       }

  // open all *R files and read basic data from them
    sprintf(fname,"%s.R",argv[1]);
    R = fopen(fname,"w"); // outpout R file
    if (R == NULL)
       { fprintf (stderr, "Cannot open %s\n",fname);
   exit(1);
       }

    // input file dicz.int.R
    sprintf(fname,"%s.dicz.int.R",argv[1]);
    diczR = fopen(fname,"r");
    if (diczR == NULL){ 
      fprintf (stderr, "Cannot open %s\n",fname); exit(1);
    }
    // get largest terminal in dicz.int.R
    largest_symbol= maxbelowUnique(diczR);
    rewind(diczR);
    // get number of terms in dicz.int.R
    e = fread(&terms,sizeof(int),1,diczR);
    die(e!=1,"Read error");
    assert(terms>Unique);  // Unique is a terminal in dicz.int
    stat (fname,&s);
    rules = (s.st_size-sizeof(int))/(2*sizeof(int));

    // input file parse.R
    sprintf(fname,"%s.parse.R",argv[1]);
    parseR = fopen(fname,"r");
    if (parseR == NULL){ 
      fprintf (stderr, "Cannot open %s\n",fname);
      exit(1);
    }
    e = fread (&phrases,sizeof(int),1,parseR);
    die(e!=1,"Read error");
    stat (fname,&s);
    prules = (s.st_size-sizeof(int))/(2*sizeof(int));

    // input file dicz.int.C copied in Adicz[]
    sprintf(fname,"%s.dicz.int.C",argv[1]);
    diczC = fopen(fname,"r");
    if (diczC == NULL){ 
      fprintf (stderr, "Cannot open %s\n",fname); exit(1);
    }
    stat(fname,&s);
    sizeC = s.st_size/sizeof(int);
    AdiczC = (int*) malloc (s.st_size);
    if (AdiczC == NULL) { 
      fprintf (stderr, "Cannot allocate %li bytes for %s.dicz.int.C\n",s.st_size,argv[1]); exit(1);
    }
    e = fread (AdiczC,sizeof(int),sizeC,diczC);
    die(e!=sizeC,"Read error");
    fclose (diczC);

    // complete computation of the largest symbol 
    for(int i=0;i<sizeC;i++)
      if(AdiczC[i]<Unique && AdiczC[i]>largest_symbol)
         largest_symbol = AdiczC[i];
    // number or terminals (upper bound), id of the first non terminal 
    v256 = largest_symbol +1;
    fprintf(stderr,"Size of the original input alphabet: %d\n",v256);

    // allocate translation table  
    transl = (int*) malloc ((phrases+1)*sizeof(int));
    if ((transl == NULL) && (phrases > 0)) { 
      fprintf (stderr, "Cannot allocate %li bytes for the phrases of %s.dicz.int.R\n", phrases*sizeof(int),argv[1]); exit(1);
    }


    // start creating .R file    
    // non terminals in the output grammar have id >=v256
    fwrite(&v256,sizeof(int),1,R); 

    // copy rules from dicz.int.R to .R shifting their index of (terms - v256)
    // leaves cursor at the end of R, ready to add more rules
    for (i=0;i<rules;i++) { 
      e = fread(val,sizeof(int),2,diczR);
      die(e!=2,"Read error");
      // right hand side of rules are true terminal or nonterminal (ie not separator)
      assert(val[0]<v256 || val[0]>=terms);
      if (val[0] >= terms) val[0] -= (terms-v256); 
      assert(val[1]<v256 || val[1]>=terms);
      if (val[1] >= terms) val[1] -= (terms-v256); 
      fwrite(val,sizeof(int),2,R);
    }
    fclose(diczR);
    assert(ftell(R)==(1+2*rules)*sizeof(int));  // written exactly "rules" rules

    // scan dicz.int.C, identify the terms and create balanced rules on them
    // adding those to R and remembering the mapping phrase -> top nonterm in transl
    i = 0; p = 1;
    while (i < sizeC) { 
      int j = i;
      while ((AdiczC[j] < v256) || (AdiczC[j] >= terms)) // until we reach a separator 
      { if (AdiczC[j] >= terms) AdiczC[j] -= (terms-v256); 
        j++;
      }
      assert(AdiczC[j]>=Unique); // this must be a separator 
      int ni = j+1; // next i 
      while (j-i > 1) { // balanced grammar construction 
        int k = i;
        int ko = i;
        while (k+1 < j) { 
          val[0] = AdiczC[k++]; val[1] = AdiczC[k++]; // new rule to generate AdiczC[k,k+1]
          fwrite(val,sizeof(int),2,R);
          AdiczC[ko++] = v256 + rules++;           // left hand side of the new rule 
        }
        if (k < j) AdiczC[ko++] = AdiczC[k];         // odd symbol  
        j = ko;
      }
      transl[p++] = AdiczC[i];  // save non terminal generating the phrase 
      i = ni;
    }
    free (AdiczC);
    assert(ftell(R)==(1+2*rules)*sizeof(int));

    // translate the rules in .parse.R according to table transl
    // terminal are transformed into their balanced grammar 
    // non terminal are shifted by  -phrases + (rules+v256)
    for (i=0;i<prules;i++)
       { e = fread(val,sizeof(int),2,parseR);
         die(e!=2,"Read error");
         if (val[0] < phrases) val[0] = transl[val[0]];   // terminal in the parse 
         else val[0] = val[0] - phrases + rules + v256;
         if (val[1] < phrases) val[1] = transl[val[1]];
         else val[1] = val[1] - phrases + rules + v256;
         fwrite(val,sizeof(int),2,R);
       }
    fclose (parseR);
    assert(ftell(R)==(1+2*(rules+prules))*sizeof(int));
    
    // output .R file completed
    if (fclose(R) != 0) { 
      fprintf (stderr,"Cannot close %s.R\n",argv[1]); exit(1);
    }
    
    // creation of the final .C file
    // translate the rules in .parse.C according to table transl
    sprintf(fname,"%s.parse.C",argv[1]);
    parseC = fopen(fname,"r");
    if (parseC == NULL) { 
      fprintf (stderr, "Cannot open %s\n",fname); exit(1);
    }
    stat (fname,&s);
    psizeC = s.st_size/sizeof(int);

    sprintf(fname,"%s.C",argv[1]);
    C = fopen(fname,"w");
    if (C == NULL) { 
      fprintf (stderr, "Cannot open %s\n",fname); exit(1);
    }

    for (i=0;i<psizeC;i++)
       { e = fread(val,sizeof(int),1,parseC);
         die(e!=1,"Read error");
         if (val[0] < phrases) val[0] = transl[val[0]];   // terminal in parse replace with grammar
         else val[0] = val[0] - phrases + rules + v256; // non terminal shifted   
         fwrite(val,sizeof(int),1,C);
       }
    free (transl);
    fclose (parseC);
    if (fclose(C) != 0) { 
      fprintf (stderr,"Cannot close %s.C\n",argv[1]); exit(1);
    }

    // estimate compression 
    if(stat(argv[1],&s)!=0) {
      fprintf(stderr,"Cannot stat file %s\n",argv[1]);
      exit(1);
    }
      
    rules += prules; // total number of rules
    long est_size = (long) ( (2.0*rules+ ((double)bits(v256+rules))*(rules+psizeC)) /8 ) +1;
    fprintf(stderr,"Prefix-Free + Repair succeeded\n");
    fprintf(stderr,"  Original file size: %li (bytes)\n",s.st_size);
    fprintf(stderr,"  Number of rules: %i\n",rules);
    fprintf(stderr,"  Final sequence length: %i\n",psizeC);
    fprintf(stderr,"  Estimated output size (bytes): %ld\n",est_size);
    fprintf(stderr,"  Compression ratio: %0.2f%%\n", (100.0* est_size)/s.st_size);
    fprintf(stdout,"  Estimated output size (stdout): %ld\n",est_size); // don't change this: est_size must be the the last printed item to stdout
    fprintf(stderr,"=== postprocessing completed!\n");
    return 0;
  }
