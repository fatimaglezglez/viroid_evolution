

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/mfe.h>
#include <math.h>
#include <stdbool.h>
#include <ViennaRNA/treedist.h>
#include <ViennaRNA/utils.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/RNAstruct.h>
#include <ViennaRNA/stringdist.h>


// NO HAY MOTIFS EN FITNESS - LEER
// qué son los unclassified de table 1?
// USE FREE()
// compute haripnins - h

double rand_uni();
long factorial(int n);
double poisson (double lambda, int k);

int main ()
{
  printf(">temporal_data: h, p, gc, L\n");
  /* definition of storage arrays */
  int N = 20; // Number of sequences // CHANGE HERE  // más de 100 me da segmentation fault
  int L = 500;
  char Nseq[N][L];
  int Nlen[N];
  char Nstructure[N][L];
  double Nmfe[N];
  int Npairs[N];
  float Ngc[N];
  int Nhairpins[N];
  double Nfitness[N];

  /* definition of selective pressures and conditions */
  float a = 2; // alfa: selección por energía (2)
  float b = 2; // beta: selección por pares de bases (1)
  float mfe_0 = -4.333; // reference energy value
  double mutations = 3.3e-3;
  double insertions = 3.3e-3;
  double deletions = 3.3e-3;

  /*** creation of N random RNA sequences ***/
  const char bases[4] = {'A', 'C', 'G', 'U'};
  char seqRNA[50];
  int n;
  for (int n=0; n<N; n++)
  {
      //srand (n*33); // SEED VALUE //  QUITAR CUANDO ESTÉ LISTO
      for (int i=0; i<30; i++)
      {
          int ri = rand_uni()*4;
          char random_base = bases[ri];
          seqRNA[i] = random_base;
      }
      strcpy(Nseq[n],seqRNA);
    }
  /*** calculation of parameters ***/
  /* CICLE */
  int T = 2; // total number of generations
  for (int t=0; t<T; t++)
  {
      //printf("\n****************************************************************************** T = %u\n", t);
      for (n=0; n<N; n++)
      {
          //printf("________________________________________________________ n = %u\n\n", n);
          /* lenght */
          int len = 0;
          for (int l=0; l<L; l++)
          {
            if ((Nseq[n][l] == 'A') || (Nseq[n][l] == 'C') || (Nseq[n][l] == 'G') || (Nseq[n][l] == 'U'))
            {len += 1;}
          }
          Nlen[n] = len;
          /* structure */
          char *structure = (char *)vrna_alloc(sizeof(char) * (sizeof(Nseq[n]) / sizeof(Nseq[n][0]) + 1));       // ?
          /* minimun folding energy */
          double mfe = vrna_circfold(Nseq[n], structure);
          strncpy(Nstructure[n], structure, sizeof(Nstructure[n]) - 1);
          Nmfe[n] = mfe;
          /* pairs */  /* G-C content */
          int p = 0;
          int pGC = 0;

          for (int l=0; l<Nlen[n]; l++)
          {
            if (structure[l] == '(')
              {p += 1;}
            if (Nseq[n][l] == 'C' || structure[l] == 'G')
              {pGC += 1;}
          }
          Npairs[n] = p;

          float ppGC = (float) pGC / (float) len;
          Ngc[n] = ppGC;
          /* fitness */
          double fit = exp(-a*mfe * (1 - mfe/(2*mfe_0)) + b*p);
          Nfitness[n] = fit;
          /*** replication ***/ // revisar para no repetir en secuencias iguales o de la misma clase
         /*printf("Sequence:\t %s\nLength:\t\t %u\nStructure:\t %s\nMFE:\t\t%6.2f\nPairs:\t\t %u\nG-C content:\t %.3f\nFitness:\t %f\n",
               Nseq[n], Nlen[n], structure, mfe, p, ppGC, fit);*/
          }
          /*hairpins*/
          for (n=0; n<N; n++)
          {
            int h = 0;
            char const* target1 = "(..";
            char const* target2 = ".)";
            const char *str = Nstructure[n];
            const char *result = str;
            while((result = strstr(result, target1)) != NULL) {
              if ((result = strstr(result, target2)) != NULL)
              {
                h += 1;
              }
              ++result;}
            printf("%s \t h: %d\n", Nstructure[n], h);
            /*char const* targetp1 = ".";
            char const* targetp2 = "..";
            char const* targetp3 = "...";
            char* pos = strstr(result, targetp3);
            if (pos) {
              if (pos - str == 0) {h += 1;}} // ...~~~~
            // FAAAAALTAAAAN CASOOOSSSSSSSSSSSSSSSSSSSSSS:  ..~~~. / .~~~.. / ~~~...*/
           Nhairpins[n] = h;


              /*{  // empieza por .
                char* pos3 = strstr(result, targetp3);
                char* pos2 = strstr(result, targetp2);
                if (pos3 - str == 0) {h += 1;} // empieza por ...
                else if (pos2 - str == Nlen[n] - 2) {h += 1;} // termina por ..
                else if ((pos2 - str == 0) && ()) {h += 1;} // termina por ..
              pos = strstr(result, targetp3);
              else if(pos - str == Nlen[n] - 3) {h += 1;}
                }
              }*/
          printf("%s \t h: %d\n", Nstructure[n], h);
          }



      double P_replication[N];
      double P_distributed[N];
      double sum_fit = 0;
      for (n=0; n<N; n++) // this is not optimal!
      {sum_fit += Nfitness[n];}
      for (n=0; n<N; n++)
      {
        P_replication[n] = Nfitness[n] / sum_fit;
        //printf("P: %f\n", P_replication[n]);
        if (n==0) {P_distributed[n] = P_replication[n];}
        else {P_distributed[n] = P_replication[n] + P_distributed[n-1];}
        //printf("D: %f\n", P_distributed[n-1]);
      }
      char Nseq_2[N][L];
      for (n=0; n<N; n++)
      {Nseq_2[n][0] = '0';}
      //srand(t*44); // REVISAR
      double r;
      int count[n];
      for (n=0; n<N; n++)
      {count[n]=0;}

      for (int rr=0; rr<N; rr++)
      {
        r = rand_uni();
        for (n=0; n<N; n++) // optimizar // pensar en agrupar en clases (misma seq)
        {
          if ((n==0 || P_distributed[n-1] <= r) && (r < P_distributed[n])) // excluyo r = 1 :(
          {
            strncpy(Nseq_2[n], Nseq[n], sizeof(Nseq_2[n]) - 1);
            count[n] += 1;
          }
        }
      }
      /*for (n=0; n<N; n++)
      {printf("Sequence:\t %s  \t\t%u\n", Nseq[n], count[n]);}*/

      int index = 0;
      for (n=0; n<N; n++)
      {
        if (Nseq_2[n][0] != '0')
        {
          for (int c=0; c<count[n]; c++)
          {
            strncpy(Nseq[index], Nseq_2[n], sizeof(Nseq[index]) -1 );
            index += 1;
          }
        }
      }

      for (n=0; n<N; n++)
      {
        /*mutations*/
        float lambda_m = mutations*Nlen[n];
        //printf ("mutations: lambda = %f \t", lambda_m);
        double r1 = rand_uni();
        int k1=-1;
        double P_poisson_sum_m[10];
        do
        {
          k1+=1;
          double pp = poisson(lambda_m, k1);
          if (k1==0) {P_poisson_sum_m[k1] = pp;}
          else {P_poisson_sum_m[k1] = pp + P_poisson_sum_m[k1-1];}
        }
        while (P_poisson_sum_m[k1] < r1);
        //printf ("Poisson sum = %lf \t k = %u\t Random number = %lf\t Sequence:\t %s\n", P_poisson_sum_m[k1], k1, r1, Nseq[n]);
        if (k1 != 0)
        {
          for (int k=0; k<k1; k++)
          {
            int r_mut = rand_uni()*Nlen[n]; // posición para la mutación
            _Bool different = true;
            while (different)
            {
              int ri1 = rand_uni()*4; // posicion en "bases" para sustituir
              char mut = bases[ri1];
              if (Nseq[n][r_mut] != mut)
              {Nseq[n][r_mut] = mut; different = false;}
            }
          }
        }
        // printf("Sequence:\t %s\n", Nseq[n]); //    ESTÁ OK SO FAR
        /*insertions*/
        float lambda_i = insertions*Nlen[n];
        double r2 = rand_uni();
        int k2=-1;
        double P_poisson_sum_i[10];
        do
        {
          k2+=1;
          double pp = poisson(lambda_i, k2);
          if (k2==0) {P_poisson_sum_i[k2] = pp;}
          else {P_poisson_sum_i[k2] = pp + P_poisson_sum_i[k2-1];}
        }
        while (P_poisson_sum_i[k2] < r2);
        //printf("k = %u\n", k2); // da 66 en mogollón de casos y en los que da 1 introduce un ">"
        //printf("Sequence:\t %s\n", Nseq[n]); // HASTA AQUÍ OK en seq
        if (k2 != 0)
        {
          for (int k=0; k<k2; k++)
          {
            int r_ins = rand_uni()*Nlen[n]; // posición para la inserción
            char temp_seq[L];
            for (int l=0; l<r_ins; l++)
            {temp_seq[l] = Nseq[n][l];}
            int ri2 = rand_uni()*4;
            //printf ("%u\t%c\n", ri2, bases[ri2]);
            temp_seq[r_ins] = bases[ri2];
            for (int l=r_ins+1; l<=L; l++) // dejamos que llegue a L porque tiene que haber un nt más
            {temp_seq[l] = Nseq[n][l-1];}
            strncpy(Nseq[n], temp_seq, sizeof(Nseq[n]) - 1);
          }
        }
        /*deletions*/
        float lambda_d = deletions*Nlen[n];
        double r3 = rand_uni();
        int k3=-1;
        double P_poisson_sum_d[10];
        do
        {
          k3+=1;
          double pp = poisson(lambda_d, k3);
          if (k3==0) {P_poisson_sum_d[k3] = pp;}
          else {P_poisson_sum_d[k3] = pp + P_poisson_sum_d[k3-1];}
        }
        while (P_poisson_sum_d[k3] < r3);
        //printf("k = %u\n", k2); // da 66 en mogollón de casos y en los que da 1 introduce un ">"
        //printf("Sequence:\t %s\n", Nseq[n]); // HASTA AQUÍ OK en seq
        if (k3 != 0)
        {
          for (int k=0; k<k3; k++)
          {
            int r_del = rand_uni()*Nlen[n]; // posición para la inserción
            char temp_seq[L];
            for (int l=0; l<r_del; l++)
            {temp_seq[l] = Nseq[n][l];}
            for (int l=r_del; l<L; l++) // dejamos que llegue a L porque tiene que haber un nt más
            {temp_seq[l] = Nseq[n][l+1];}
            strncpy(Nseq[n], temp_seq, sizeof(Nseq[n]) - 1);
          }
        }
      //printf("\t %s\n", Nseq[n]);
    } // cierra n


    double sum_h = 0;
    double sum_p = 0;
    double sum_gc = 0;
    double sum_L = 0;
    double mean_h;
    double mean_p;
    double mean_gc;
    double mean_L;
    for (int n=0; n<N; n++)
    {
      sum_h += Nhairpins[n];
      sum_p += Npairs[n];
      sum_gc += Ngc[n];
      sum_L += Nlen[n];
    }
    mean_h = sum_h / (double) N;
    mean_p = sum_p / (double) N;
    mean_gc = sum_gc / (double) N;
    mean_L = sum_L / (double) N;

    printf("[%f, %f, %f, %f],", mean_h, mean_p, mean_gc, mean_L);
  }    // cierra el ciclo t


  printf("\n");
  printf(">N\n%u\n", N);
  printf(">T\n%u\n", T);
  printf(">a\n%.2f\n", a);
  printf(">b\n%.2f", b);
  printf("\n>sequences\n");
  for (int n=0; n<N; n++)
  {printf("%s,", Nseq[n]);}
  printf("\n>lenghts\n");
  for (int n=0; n<N; n++)
  {printf("%u,", Nlen[n]);}
  printf("\n>structures\n");
  for (int n=0; n<N; n++)
  {printf("%s,", Nstructure[n]);}
  printf("\n>MFE\n");
  for (int n=0; n<N; n++)
  {printf("%f,", Nmfe[n]);}
  printf("\n>pairs\n");
  for (int n=0; n<N; n++)
  {printf("%u,", Npairs[n]);}
  printf("\n>GC_content\n");
  for (int n=0; n<N; n++)
  {printf("%f,", Ngc[n]);}
  printf("\n");

  return 0;
} // cierra main

double rand_uni() {return (double) rand() / (double) RAND_MAX;}

long factorial(int n)
{
  int c;
  long f = 1;
  for (c = 1; c <= n; c++)
    {f = f * c;}
  return f;
}

double poisson(double lambda, int k) {return (pow(lambda, k) * exp(-lambda)) / factorial(k);}
