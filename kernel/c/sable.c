#include "easypap.h"

#include <omp.h>
#include <stdbool.h>
#include <sys/mman.h>
#include <unistd.h>

typedef unsigned TYPE;

static TYPE *TABLE = NULL;
static TYPE *STABILITY_TABLE = NULL;

static volatile int changement;

static TYPE max_grains;

static inline TYPE *table_cell(TYPE *restrict i, int y, int x)
{
  return i + y * DIM + x;
}

static inline TYPE *stable_table_cell(TYPE *restrict i, int y, int x)
{
  return i + y * (DIM / TILE_W) + x;
}

#define TILE_SIZE TILE_W *TILE_H

#define table(y, x) (*table_cell(TABLE, (y), (x)))
#define not_stable(y, x) (*stable_table_cell(STABILITY_TABLE, (y), (x)))

#define RGB(r, g, b) rgba(r, g, b, 0xFF)

////////////////////// DEBUG ////////////////////////
void print_debug(TYPE *t, int nth, int ntw)
{
  for (int i = 0; i < nth; i++)
  {
    printf("%d : ", i);
    for (int j = 0; j < ntw; j++)
    {
      printf("%u | ", *stable_table_cell(t, i, j));
    }
    printf("\n");
  }
}

void sable_init()
{
  if (TABLE == NULL)
  {
    const unsigned size = DIM * DIM * sizeof(TYPE);
    const unsigned stability_size = (DIM / TILE_H) * (DIM / TILE_W) * sizeof(TYPE);

    PRINT_DEBUG('u', "Memory footprint = 2 x %d bytes\n", size);

    TABLE = mmap(NULL, size, PROT_READ | PROT_WRITE,
                 MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

    STABILITY_TABLE = mmap(NULL, stability_size, PROT_READ | PROT_WRITE,
                           MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
        not_stable(y / TILE_H, x / TILE_W) = 1;
  }
}
void sable_finalize()
{
  const unsigned size = DIM * DIM * sizeof(TYPE);
  const unsigned stability_size = (DIM / TILE_H) * (DIM / TILE_W) * sizeof(TYPE);

  munmap(TABLE, size);
  munmap(STABILITY_TABLE, stability_size);
}

///////////////////////////// Production d'une image
void sable_refresh_img()
{
  unsigned long int max = 0;
  for (int i = 1; i < DIM - 1; i++)
    for (int j = 1; j < DIM - 1; j++)
    {
      int g = table(i, j);
      int r, v, b;
      r = v = b = 0;
      if (g == 1)
        v = 255;
      else if (g == 2)
        b = 255;
      else if (g == 3)
        r = 255;
      else if (g == 4)
        r = v = b = 255;
      else if (g > 4)
        r = b = 255 - (240 * ((double)g) / (double)max_grains);

      cur_img(i, j) = RGB(r, v, b);
      if (g > max)
        max = g;
    }
  max_grains = max;
}

///////////////////////////// Version séquentielle simple (seq)

static inline int compute_new_state(int y, int x)
{
  if (table(y, x) >= 4)
  {
    unsigned long int div4 = table(y, x) / 4;
    table(y, x - 1) += div4;
    table(y, x + 1) += div4;
    table(y - 1, x) += div4;
    table(y + 1, x) += div4;
    table(y, x) %= 4;
    return 1;
  }
  return 0;
}

static int do_tile(int x, int y, int width, int height, int who)
{
  int chgt = 0;
  PRINT_DEBUG('c', "tuile [%d-%d][%d-%d] traitée\n", x, x + width - 1, y,
              y + height - 1);

  monitoring_start_tile(who);

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
    {
      chgt |= compute_new_state(i, j);
    }

  monitoring_end_tile(x, y, width, height, who);
  return chgt;
}
static int do_double_tile(int x, int y, int width, int height, int who)
{
  int chgt = 0;
  PRINT_DEBUG('c', "tuile [%d-%d][%d-%d] traitée\n", x, x + width - 1, y,
              y + height - 1);

  monitoring_start_tile(who);

  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
    {
      chgt |= compute_new_state(i, j);
    }
  if (chgt == 0)
  {
    monitoring_end_tile(x, y, width, height, who);
    return 0;
  }
  for (int i = y; i < y + height; i++)
    for (int j = x; j < x + width; j++)
    {
      chgt |= compute_new_state(i, j);
    }
  monitoring_end_tile(x, y, width, height, who);
  return chgt;
}
static int do_tile_stable(int x, int y, int width, int height, int who)
{
  int change = 0;
  monitoring_start_tile(who);
  for (int i = y; i < y + height; i++)
  {
    change += compute_new_state(i, x);
    change += compute_new_state(i, x + width - 1);
  }
  for (int j = x; j < x + width; j++)
  {
    change += compute_new_state(y, j);
    change += compute_new_state(y + height - 1, j);
  }
  monitoring_end_tile(x, y, width, height, who);
  if (change == 0)
  {
    return 0;
  }
  return 1;
}
static int do_tile_unstable(int x, int y, int width, int height, int who)
{
  int change = 0;
  PRINT_DEBUG('c', "tuile [%d-%d][%d-%d] traitée\n", x, x + width - 1, y,
              y + height - 1);

  monitoring_start_tile(who);
  //the  tile is unstable so we need to check all pixels
  for (int i = y; i < y + height; i++)
  {
    for (int j = x; j < x + width; j++)
    {
      change += compute_new_state(i, j);
    }
  }
  monitoring_end_tile(x, y, width, height, who);
  if (!change)
  {
    return 0;
  }
  return 1;
}
static int do_double_tile_stable(int x, int y, int width, int height, int who)
{
  int change = 0;
  monitoring_start_tile(who);
  for (int i = y; i < y + height; i++)
  {
    change += compute_new_state(i, x);
    change += compute_new_state(i, x + width - 1);
  }
  for (int j = x; j < x + width; j++)
  {
    change += compute_new_state(y, j);
    change += compute_new_state(y + height - 1, j);
  }
  if (change == 0)
  {
    monitoring_end_tile(x, y, width, height, who);
    return 0;
  }
  change += do_tile_unstable(x, y, width, height, who);
  monitoring_end_tile(x, y, width, height, who);
  if (change == 0)
    return 0;
  return 1;
}
static int do_double_tile_unstable(int x, int y, int width, int height, int who)
{
  int change = 0;
  PRINT_DEBUG('c', "tuile [%d-%d][%d-%d] traitée\n", x, x + width - 1, y,
              y + height - 1);

  monitoring_start_tile(who);
  //the  tile is unstable so we need to check all pixels
  for (int i = y; i < y + height; i++)
  {
    for (int j = x; j < x + width; j++)
    {
      change += compute_new_state(i, j);
    }
  }
  if (!change)
  {
    monitoring_end_tile(x, y, width, height, who);
    return 0;
  }

  for (int i = y; i < y + height; i++)
  {
    for (int j = x; j < x + width; j++)
    {
      change += compute_new_state(i, j);
    }
  }
  monitoring_end_tile(x, y, width, height, who);
  if (!change)
  {
    return 0;
  }
  return 1;
}
unsigned sable_compute_seq(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    changement = 0;
    // On traite toute l'image en un coup (oui, c'est une grosse tuile)
    changement |= do_tile(1, 1, DIM - 2, DIM - 2, 0);
    if (changement == 0)
      return it;
  }
  return 0;
}

///////////////////////////// Version tuilée (tiled)

unsigned sable_compute_tiled(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    changement = 0;

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
        changement |= do_tile(x + (x == 0), y + (y == 0),
                              TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                              TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                              0 /* CPU id */);
    if (changement == 0)
      return it;
  }

  return 0;
}

unsigned sable_compute_double_tiled(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    changement = 0;

    for (int y = 0; y < DIM; y += TILE_H)
      for (int x = 0; x < DIM; x += TILE_W)
        changement |= do_double_tile(x + (x == 0), y + (y == 0),
                                     TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                     TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                     0 /* CPU id */);
    if (changement == 0)
      return it;
  }

  return 0;
}

unsigned sable_compute_tiled_stable(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    changement = 0;

    for (int i = 0; i < DIM; i += TILE_H)
    {
      for (int j = 0; j < DIM; j += TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_unstable(j + (j == 0),
                                                                i + (i == 0),
                                                                TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_stable(j + (j == 0),
                                                              i + (i == 0),
                                                              TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                              TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                              omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }
    if (changement == 0)
      return it;
  }
  return 0;
}

unsigned sable_compute_double_tiled_stable(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    changement = 0;

    for (int i = 0; i < DIM; i += TILE_H)
    {
      for (int j = 0; j < DIM; j += TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_unstable(j + (j == 0),
                                                                       i + (i == 0),
                                                                       TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                       TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                       omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_stable(j + (j == 0),
                                                                     i + (i == 0),
                                                                     TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                     TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                     omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }
    if (changement == 0)
      return it;
  }
  return 0;
}

////////////////////////////// Version OMP
unsigned sable_compute_tiled_omp(unsigned nb_iter) //tile size perfect is 32
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    changement = 0;

#pragma omp parallel for collapse(2) reduction(| \
                                               : changement) schedule(runtime)
    for (int y = 0; y < DIM; y += 2 * TILE_H)
      for (int x = 0; x < DIM; x += 2 * TILE_W)
        changement |= do_tile(x + (x == 0), y + (y == 0),
                              TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                              TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                              omp_get_thread_num());

#pragma omp parallel for collapse(2) reduction(| \
                                               : changement) schedule(runtime)
    for (int y = TILE_H; y < DIM; y += 2 * TILE_H)
      for (int x = TILE_W; x < DIM; x += 2 * TILE_W)
        changement |= do_tile(x + (x == 0), y + (y == 0),
                              TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                              TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                              omp_get_thread_num());

#pragma omp parallel for collapse(2) reduction(| \
                                               : changement) schedule(runtime)
    for (int y = 0; y < DIM; y += 2 * TILE_H)
      for (int x = TILE_W; x < DIM; x += 2 * TILE_W)
        changement |= do_tile(x + (x == 0), y + (y == 0),
                              TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                              TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                              omp_get_thread_num());

#pragma omp parallel for collapse(2) reduction(| \
                                               : changement) schedule(runtime)
    for (int y = TILE_H; y < DIM; y += 2 * TILE_H)
      for (int x = 0; x < DIM; x += 2 * TILE_W)
        changement |= do_tile(x + (x == 0), y + (y == 0),
                              TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                              TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                              omp_get_thread_num());
    if (changement == 0)
      return it;
  }

  return 0;
}

unsigned sable_compute_double_tiled_omp(unsigned nb_iter) //tile size perfect is 32
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    changement = 0;

#pragma omp parallel for collapse(2) reduction(| \
                                               : changement) schedule(runtime)
    for (int y = 0; y < DIM; y += 2 * TILE_H)
      for (int x = 0; x < DIM; x += 2 * TILE_W)
        changement |= do_double_tile(x + (x == 0), y + (y == 0),
                                     TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                     TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                     omp_get_thread_num());

#pragma omp parallel for collapse(2) reduction(| \
                                               : changement) schedule(runtime)
    for (int y = TILE_H; y < DIM; y += 2 * TILE_H)
      for (int x = TILE_W; x < DIM; x += 2 * TILE_W)
        changement |= do_double_tile(x + (x == 0), y + (y == 0),
                                     TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                     TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                     omp_get_thread_num());

#pragma omp parallel for collapse(2) reduction(| \
                                               : changement) schedule(runtime)
    for (int y = 0; y < DIM; y += 2 * TILE_H)
      for (int x = TILE_W; x < DIM; x += 2 * TILE_W)
        changement |= do_double_tile(x + (x == 0), y + (y == 0),
                                     TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                     TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                     omp_get_thread_num());

#pragma omp parallel for collapse(2) reduction(| \
                                               : changement) schedule(runtime)
    for (int y = TILE_H; y < DIM; y += 2 * TILE_H)
      for (int x = 0; x < DIM; x += 2 * TILE_W)
        changement |= do_double_tile(x + (x == 0), y + (y == 0),
                                     TILE_W - ((x + TILE_W == DIM) + (x == 0)),
                                     TILE_H - ((y + TILE_H == DIM) + (y == 0)),
                                     omp_get_thread_num());
    if (changement == 0)
      return it;
  }

  return 0;
}

unsigned sable_compute_tiled_stable_omp(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    changement = 0;

#pragma omp parallel for collapse(2) reduction(+ \
                                               : changement) schedule(runtime)
    for (int i = 0; i < DIM; i += 2 * TILE_H)
    {
      for (int j = 0; j < DIM; j += 2 * TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_unstable(j + (j == 0),
                                                                i + (i == 0),
                                                                TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_stable(j + (j == 0),
                                                              i + (i == 0),
                                                              TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                              TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                              omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }

#pragma omp parallel for collapse(2) reduction(+ \
                                               : changement) schedule(runtime)
    for (int i = TILE_H; i < DIM; i += 2 * TILE_H)
    {
      for (int j = 0; j < DIM; j += 2 * TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_unstable(j + (j == 0),
                                                                i + (i == 0),
                                                                TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_stable(j + (j == 0),
                                                              i + (i == 0),
                                                              TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                              TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                              omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }

#pragma omp parallel for collapse(2) reduction(+ \
                                               : changement) schedule(runtime)
    for (int i = 0; i < DIM; i += 2 * TILE_H)
    {
      for (int j = TILE_W; j < DIM; j += 2 * TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_unstable(j + (j == 0),
                                                                i + (i == 0),
                                                                TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_stable(j + (j == 0),
                                                              i + (i == 0),
                                                              TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                              TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                              omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }

#pragma omp parallel for collapse(2) reduction(+ \
                                               : changement) schedule(runtime)
    for (int i = TILE_H; i < DIM; i += 2 * TILE_H)
    {
      for (int j = TILE_W; j < DIM; j += 2 * TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_unstable(j + (j == 0),
                                                                i + (i == 0),
                                                                TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_tile_stable(j + (j == 0),
                                                              i + (i == 0),
                                                              TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                              TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                              omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }
    if (changement == 0)
      return it;
  }
  return 0;
}

unsigned sable_compute_double_tiled_stable_omp(unsigned nb_iter)
{
  for (unsigned it = 1; it <= nb_iter; it++)
  {
    changement = 0;

#pragma omp parallel for collapse(2) reduction(+ \
                                               : changement) schedule(runtime)
    for (int i = 0; i < DIM; i += 2 * TILE_H)
    {
      for (int j = 0; j < DIM; j += 2 * TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_unstable(j + (j == 0),
                                                                       i + (i == 0),
                                                                       TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                       TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                       omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_stable(j + (j == 0),
                                                                     i + (i == 0),
                                                                     TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                     TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                     omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }

#pragma omp parallel for collapse(2) reduction(+ \
                                               : changement) schedule(runtime)
    for (int i = TILE_H; i < DIM; i += 2 * TILE_H)
    {
      for (int j = 0; j < DIM; j += 2 * TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_unstable(j + (j == 0),
                                                                       i + (i == 0),
                                                                       TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                       TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                       omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_stable(j + (j == 0),
                                                                     i + (i == 0),
                                                                     TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                     TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                     omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }

#pragma omp parallel for collapse(2) reduction(+ \
                                               : changement) schedule(runtime)
    for (int i = 0; i < DIM; i += 2 * TILE_H)
    {
      for (int j = TILE_W; j < DIM; j += 2 * TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_unstable(j + (j == 0),
                                                                       i + (i == 0),
                                                                       TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                       TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                       omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_stable(j + (j == 0),
                                                                     i + (i == 0),
                                                                     TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                     TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                     omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }

#pragma omp parallel for collapse(2) reduction(+ \
                                               : changement) schedule(runtime)
    for (int i = TILE_H; i < DIM; i += 2 * TILE_H)
    {
      for (int j = TILE_W; j < DIM; j += 2 * TILE_W)
      {
        if (not_stable(i / TILE_H, j / TILE_W) || j == 0 || i == 0 || i + TILE_H >= DIM || j + TILE_W >= DIM)
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_unstable(j + (j == 0),
                                                                       i + (i == 0),
                                                                       TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                       TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                       omp_get_thread_num());
        }
        else if (not_stable((i / TILE_H) - 1, j / TILE_W) || not_stable((i / TILE_H) + 1, j / TILE_W) || not_stable((i / TILE_H), (j / TILE_W) - 1) || not_stable((i / TILE_H), (j / TILE_W) + 1))
        {
          not_stable(i / TILE_H, j / TILE_W) = do_double_tile_stable(j + (j == 0),
                                                                     i + (i == 0),
                                                                     TILE_W - ((j + TILE_W == DIM) + (j == 0)),
                                                                     TILE_H - ((i + TILE_H == DIM) + (i == 0)),
                                                                     omp_get_thread_num());
        }
        changement += not_stable(i / TILE_H, j / TILE_W);
      }
    }
    if (changement == 0)
      return it;
  }
  return 0;
}

///////////////////////////// open CL
static cl_mem changed, ocl_changes;
static cl_mem in_state, out_state;
static volatile int changement;

void swap_mem(cl_mem *m1, cl_mem *m2)
{
  cl_mem tmp = *m1;
  *m1 = *m2;
  *m2 = tmp;
}

void sable_init_ocl(void)
{
  const unsigned size = DIM * DIM * sizeof(TYPE);

  PRINT_DEBUG('u', "Memory footprint = 2 x %d bytes\n", size);

  TABLE = mmap(NULL, size, PROT_READ | PROT_WRITE,
               MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

  changed = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, NULL);
  if (!changed)
    exit_with_error("Failed to allocate ocl changes variable");
}

void sable_refresh_img_ocl()
{
  cl_int err;
  err = clEnqueueReadBuffer(queue, cur_buffer, CL_TRUE, 0,
                            sizeof(TYPE) * DIM * DIM, TABLE, 0, NULL, NULL);
  check(err, "Failed to read buffer from GPU");

  sable_refresh_img();
}

unsigned sable_invoke_ocl(unsigned nb_iter)
{
  int chgt;
  size_t global[2] = {GPU_SIZE_X, GPU_SIZE_Y};
  size_t local[2] = {GPU_TILE_W, GPU_TILE_H};
  cl_int err;

  monitoring_start_tile(easypap_gpu_lane(TASK_TYPE_COMPUTE));

  for (unsigned it = 1; it <= nb_iter; it++)
  {
    chgt = 0;
    if (it % (refresh_rate - 1) == 0)
    {
      check(
          clEnqueueWriteBuffer(queue, changed, CL_TRUE, 0,
                               sizeof(int), &chgt, 0, NULL, NULL),
          "Failed to write in changed");
    }

    err = 0;
    err |= clSetKernelArg(compute_kernel, 0, sizeof(cl_mem), &cur_buffer);
    err |= clSetKernelArg(compute_kernel, 1, sizeof(cl_mem), &next_buffer);
    err |= clSetKernelArg(compute_kernel, 2, sizeof(cl_mem), &changed);
    check(err, "Failed to set kernel arguments");

    err = clEnqueueNDRangeKernel(queue, compute_kernel, 2, NULL, global, local,
                                 0, NULL, NULL);
    check(err, "Failed to execute kernel");
    if (it % (refresh_rate - 1) == 0)
    {
      clFinish(queue);
      check(
          clEnqueueReadBuffer(queue, changed, CL_TRUE, 0,
                              sizeof(int), &chgt, 0, NULL, NULL),
          "Failed to read in changed");
    }

    {
      if (it % (refresh_rate - 1) == 0 && chgt == 0)
      {
        return it;
      }
      swap_mem(&cur_buffer, &next_buffer);
    }
  }
  monitoring_end_tile(0, 0, DIM, DIM, easypap_gpu_lane(TASK_TYPE_COMPUTE));

  return 0;
}
void sable_init_ocl_freq(void)
{
  const unsigned size = DIM * DIM * sizeof(TYPE);

  TABLE = mmap(NULL, size, PROT_READ | PROT_WRITE,
               MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

  ocl_changes = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), NULL, NULL);
  if (!ocl_changes)
    exit_with_error("Failed to allocate cl_mem variable: ocl_changes");
}
void sable_refresh_img_ocl_freq()
{
  cl_int err = clEnqueueReadBuffer(queue, cur_buffer, CL_TRUE, 0,
                                   sizeof(TYPE) * DIM * DIM, TABLE, 0, NULL, NULL);
  check(err, "Failed to read buffer from GPU");
  sable_refresh_img();
}

unsigned sable_invoke_ocl_freq(unsigned nb_iter)
{
  // local domain size for our calculation
  size_t global[2] = {GPU_SIZE_X, GPU_SIZE_Y};
  size_t local[2] = {GPU_TILE_W, GPU_TILE_H};
  cl_int err;
  int changement;
  int frequence_check_changes = refresh_rate;

  monitoring_start_tile(easypap_gpu_lane(TASK_TYPE_COMPUTE));

  for (unsigned it = 1; it <= nb_iter; it++)
  {
    if (it % frequence_check_changes == 0)
    {
      // on écrit 0 dans la détection des changements
      changement = 0;
      check(
          clEnqueueWriteBuffer(queue, ocl_changes, CL_TRUE, 0,
                               sizeof(int), &changement, 0, NULL, NULL),
          "Failed to write in  ocl_changes");
    }
    // Set kernel arguments
    //
    err = 0;
    err |= clSetKernelArg(compute_kernel, 0, sizeof(cl_mem), &cur_buffer);
    err |= clSetKernelArg(compute_kernel, 1, sizeof(cl_mem), &next_buffer);
    err |= clSetKernelArg(compute_kernel, 2, sizeof(cl_mem), &ocl_changes);
    check(err, "Failed to set kernel arguments");

    err = clEnqueueNDRangeKernel(queue, compute_kernel, 2, NULL, global, local,
                                 0, NULL, NULL);
    check(err, "Failed to execute kernel");
    if (it % frequence_check_changes == 0)
    {
      clFinish(queue);
      check(
          clEnqueueReadBuffer(queue, ocl_changes, CL_TRUE, 0,
                              sizeof(int), &changement, 0, NULL, NULL),
          "Failed to read in changement");
      if (changement == 0)
      {
        return it;
      }
      swap_mem(&cur_buffer, &next_buffer);
    }
  }
  monitoring_end_tile(0, 0, DIM, DIM, easypap_gpu_lane(TASK_TYPE_COMPUTE));
  return 0;
}

///////////////////////////// Configurations initiales

static void sable_draw_4partout(void);

void sable_draw(char *param)
{
  // Call function ${kernel}_draw_${param}, or default function (second
  // parameter) if symbol not found
  hooks_draw_helper(param, sable_draw_4partout);
}

void sable_draw_4partout(void)
{
  max_grains = 8;
  for (int i = 1; i < DIM - 1; i++)
    for (int j = 1; j < DIM - 1; j++)
      cur_img(i, j) = table(i, j) = 4;
}

void sable_draw_DIM(void)
{
  max_grains = DIM;
  for (int i = DIM / 4; i < DIM - 1; i += DIM / 4)
    for (int j = DIM / 4; j < DIM - 1; j += DIM / 4)
      cur_img(i, j) = table(i, j) = i * j / 4;
}

void sable_draw_alea(void)
{
  max_grains = 5000;
  for (int i = 0; i < DIM >> 3; i++)
  {
    int i = 1 + random() % (DIM - 2);
    int j = 1 + random() % (DIM - 2);
    int grains = 1000 + (random() % (4000));
    cur_img(i, j) = table(i, j) = grains;
  }
}