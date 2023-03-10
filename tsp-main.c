#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>


// #include <omp.h>

#define DEBUG_PROGRESS

#define MAX_NBVILLES 22

typedef int DTab_t[MAX_NBVILLES][MAX_NBVILLES];
typedef int chemin_t[MAX_NBVILLES];

/* macro de mesure de temps, retourne une valeur en �secondes */
#define TIME_DIFF(t1, t2) \
  ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

/* dernier minimum trouv� */
int minimum = INT_MAX;

/* tableau des distances */
DTab_t distance;

/* nombre de villes */
int nbVilles;

/* profondeur du parallélisme */
int grain;

#define MAXX 100
#define MAXY 100
typedef struct
{
  int x, y;
} coor_t;

typedef coor_t coortab_t[MAX_NBVILLES];

#include "collapse.c"

void initialisation(int Argc, char *Argv[])
{

  if (Argc < 4 || Argc > 5)
  {
    fprintf(stderr, "Usage: %s  <nbVilles> <seed> [grain] <kernel>\n", Argv[0]);
    exit(1);
  }

 grain = (Argc == 5) ? atoi(Argv[3]) : 0;


  /* initialisation du tableau des distances */
  /* on positionne les villes aléatoirement sur une carte MAXX x MAXY  */

  coortab_t lesVilles;

  int i, j;
  int dx, dy;

  nbVilles = atoi(Argv[1]);
  if (nbVilles > MAX_NBVILLES)
  {
    fprintf(stderr, "trop de villes, augmentez MAX_NBVILLES\n");
    exit(1);
  }

  srand(atoi(Argv[2]));

  for (i = 0; i < nbVilles; i++)
  {
    lesVilles[i].x = rand() % MAXX;
    lesVilles[i].y = rand() % MAXY;
  }

  for (i = 0; i < nbVilles; i++)
    for (j = 0; j < nbVilles; j++)
    {
      dx = lesVilles[i].x - lesVilles[j].x;
      dy = lesVilles[i].y - lesVilles[j].y;
      distance[i][j] = (int)sqrt((double)((dx * dx) + (dy * dy)));
    }
}

/* résolution du problème du voyageur de commerce */

inline int present(int ville, int mask)
{
  return mask & (1 << ville);
}

void verifier_minimum(int lg, chemin_t chemin)
{
  int d = lg + distance[0][chemin[nbVilles - 1]];
  if (d < minimum) {
    #pragma omp critical
    {
      if (d < minimum)
      {
        minimum = d;
#ifdef DEBUG_PROGRESS
        printf("%3d :", minimum);
        for (int i = 0; i < nbVilles; i++)
          printf("%2d ", chemin[i]);
        printf("\n");
#endif
      }
    }
  }
}

void tsp_seq(int etape, int lg, chemin_t chemin, int mask)
{
  int ici, dist;

  if (etape == nbVilles)
    verifier_minimum(lg, chemin);
  else
  {
    ici = chemin[etape - 1];

    for (int i = 1; i < nbVilles; i++)
    {
      if (!present(i, mask))
      {
        chemin[etape] = i;
        dist = distance[ici][i];
        tsp_seq(etape + 1, lg + dist, chemin, mask | (1 << i));
      }
    }
  }
}

void tsp_ompfor(int etape, int lg, chemin_t chemin, int mask)
{

  if (etape > grain) { // version séquentielle
    tsp_seq(etape, lg, chemin, mask);}
    
  else { // version parallèle

    if (etape == nbVilles)
      verifier_minimum(lg, chemin);

    int ici, dist;
    ici = chemin[etape - 1];

    #pragma omp parallel firstprivate(dist)
    {
      chemin_t monChemin;
      memcpy(monChemin, chemin, sizeof(chemin_t));
      #pragma omp for 
      for (int i = 1; i < nbVilles; i++)
      {
        if (!present(i, mask))
        {
          monChemin[etape] = i;
          dist = distance[ici][i];
          tsp_ompfor(etape + 1, lg + dist, monChemin, mask | (1 << i));
        }
      }
    }
  }
}

void tsp_seq_optimized(int etape, int lg, chemin_t chemin, int mask)
{
  if (lg + distance[0][chemin[etape-1]]>= minimum)
    return;

  int ici, dist;

  if (etape == nbVilles)
    verifier_minimum(lg, chemin);
  else
  {
    ici = chemin[etape - 1];

    for (int i = 1; i < nbVilles; i++)
    {
      if (!present(i, mask))
      {
        chemin[etape] = i;
        dist = distance[ici][i];
        tsp_seq_optimized(etape + 1, lg + dist, chemin, mask | (1 << i));
      }
    }
  }
}

void tsp_ompfor_optimized(int etape, int lg, chemin_t chemin, int mask)
{

  if (etape > grain) { // version séquentielle
    tsp_seq_optimized(etape, lg, chemin, mask);}
    
  else { // version parallèle

    if (etape == nbVilles)
      verifier_minimum(lg, chemin);

    int ici, dist;
    ici = chemin[etape - 1];

    #pragma omp parallel firstprivate(dist)
    {
      chemin_t monChemin;
      memcpy(monChemin, chemin, sizeof(chemin_t));

      #pragma omp for 
      for (int i = 1; i < nbVilles; i++)
      {
        if (!present(i, mask))
        {
          monChemin[etape] = i;
          dist = distance[ici][i];
          tsp_ompfor_optimized(etape + 1, lg + dist, monChemin, mask | (1 << i));
        }
      }
    }
  }
}

void tsp_task_optimized(int etape, int lg, chemin_t chemin, int mask)
{

  if (etape > grain) { // version séquentielle
    tsp_seq_optimized(etape, lg, chemin, mask);}
    
  else { // version parallèle

    if (etape == nbVilles)
      verifier_minimum(lg, chemin);

    int ici, dist;
    ici = chemin[etape - 1];

    for (int i = 1; i < nbVilles; i++)
    {
      int *monChemin = malloc(nbVilles * sizeof(int));
      if (!monChemin) {
        fprintf(stderr, "malloc\n");
      }
      memcpy(monChemin, chemin, nbVilles * sizeof(int));

      #pragma omp task firstprivate(ici, i, mask) shared(monChemin)
      {
        if (!present(i, mask))
        {
          monChemin[etape] = i;
          dist = distance[ici][i];
          tsp_task_optimized(etape + 1, lg + dist, monChemin, mask | (1 << i));
        }
        free(monChemin);
      }

      
    }
  }
}

int main(int argc, char **argv)
{
  unsigned long temps;
  struct timeval t1, t2;

  initialisation(argc, argv);

  printf("nbVilles = %3d - grain %d \n", nbVilles, grain);
  omp_set_nested(1);

  if (!strcmp(argv[argc-1],"compare")) {
    for (int i = 0; i < 5; ++i) {
      chemin_t chemin;
      minimum = INT_MAX;
      chemin[0] = 0;
      char *fct;
      gettimeofday(&t1, NULL);
      switch (i) {
        case 0: fct = "ompcol4_opt"; printf("======= %s =======\n", fct); tsp_ompcol4_opt(); break;
        case 1: fct = "ompcol3_opt"; printf("======= %s =======\n", fct); tsp_ompcol3_opt(); break;
        case 2: fct = "ompcol2_opt"; printf("======= %s =======\n", fct); tsp_ompcol2_opt(); break;
        case 3: fct = "seq_opt"; printf("======= %s =======\n", fct); tsp_seq_optimized(1, 0, chemin, 1); break;
        case 4: fct = "ompfor_opt"; printf("======= %s =======\n", fct); tsp_ompfor_optimized(1, 0, chemin, 1); break;
      }
      gettimeofday(&t2, NULL);
      temps = TIME_DIFF(t1, t2);
      fprintf(stderr, "%s: %ld.%03ld\n", fct, temps / 1000, temps % 1000);
    }

  } else {
    chemin_t chemin;
    //omp_set_max_active_levels(grain);
  
    gettimeofday(&t1, NULL);

    chemin[0] = 0;

    if (!strcmp(argv[argc-1],"seq")) {
      tsp_seq(1, 0, chemin, 1);
    } else if (!strcmp(argv[argc-1],"seq_opt")) {
      tsp_seq_optimized(1, 0, chemin, 1);
    } else if (!strcmp(argv[argc-1],"ompfor")) {
      tsp_ompfor(1, 0, chemin, 1);
    } else if (!strcmp(argv[argc-1],"ompfor_opt")) {
      tsp_ompfor_optimized(1, 0, chemin, 1);
    } else if (!strcmp(argv[argc-1],"ompcol2")) {
      tsp_ompcol2_opt(1, 0, chemin, 1);
    } else if (!strcmp(argv[argc-1],"ompcol3")) {
      tsp_ompcol3_opt(1, 0, chemin, 1);
    } else if (!strcmp(argv[argc-1],"ompcol4")) {
      tsp_ompcol4_opt(1, 0, chemin, 1);
    } else if (!strcmp(argv[argc-1],"task_opt")){ 
      #pragma omp parallel master
      tsp_task_optimized(1, 0, chemin, 1);
    }
    else
    {
      printf("kernel inconnu\n");
      exit(1);
    }
      
    gettimeofday(&t2, NULL);

    temps = TIME_DIFF(t1, t2);
    fprintf(stderr, "%ld.%03ld\n", temps / 1000, temps % 1000);
  }
  return 0;
}
