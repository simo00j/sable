#include "kernel/ocl/common.cl"

// DO NOT MODIFY: this kernel updates the OpenGL texture buffer
// This is a sable-specific version (generic version is defined in common.cl)
__kernel void sable_update_texture(__global unsigned *cur,
                                   __write_only image2d_t tex) {
  int y = get_global_id(1);
  int x = get_global_id(0);
  int2 pos = (int2)(x, y);
  unsigned c = cur[y * DIM + x];
  unsigned r = 0, v = 0, b = 0;

  if (c == 1)
    v = 255;
  else if (c == 2)
    b = 255;
  else if (c == 3)
    r = 255;
  else if (c == 4)
    r = v = b = 255;
  else if (c > 4)
    r = v = b = (2 * c);

  c = rgba(r, v, b, 0xFF);

  write_imagef(tex, pos, color_scatter(c));
}

__kernel void sable_ocl(__global unsigned *in, __global unsigned *out,
                        __global int *changed) {
  int pix = 0;
  const int x = get_global_id(0);
  const int y = get_global_id(1);

  pix = in[y * GPU_SIZE_X + x] % 4;

  if (x + 1 < GPU_SIZE_X - 1)
    pix += in[y * GPU_SIZE_X + (x + 1)] >> 2;
  if (x - 1 > 0)
    pix += in[y * GPU_SIZE_X + (x - 1)] >> 2;

  if (y + 1 < GPU_SIZE_Y - 1)
    pix += in[(y + 1) * GPU_SIZE_X + x] >> 2;
  if (y - 1 > 0)
    pix += in[(y - 1) * GPU_SIZE_X + x] >> 2;

  barrier(CLK_LOCAL_MEM_FENCE);

  if (pix != in[y * GPU_SIZE_X + x]) {
    *changed = 1;
  }
  out[y * GPU_SIZE_X + x] = pix;
}

__kernel void sable_ocl_freq(__global unsigned *in, __global unsigned *out,
                             __global int *changes) {

  const int x = get_global_id(0);
  const int y = get_global_id(1);

  int tmp = in[y * GPU_SIZE_X + x] % 4;

  // Update gauche/droite
  tmp += (x + 1 < GPU_SIZE_X - 1) ? in[y * GPU_SIZE_X + (x + 1)] >> 2 : 0;
  tmp += (x - 1 > 0) ? in[y * GPU_SIZE_X + (x - 1)] >> 2 : 0;
  tmp += (y + 1 < GPU_SIZE_Y - 1) ? in[(y + 1) * GPU_SIZE_X + x] >> 2 : 0;
  tmp += (y - 1 > 0) ? in[(y - 1) * GPU_SIZE_X + x] >> 2 : 0;

  // barrier(CLK_LOCAL_MEM_FENCE);

  // Si on a fait un changement on met changes  1
  if (*changes == 0 && out[y * DIM + x] != in[y * DIM + x])
    *changes = 1;
  out[y * GPU_SIZE_X + x] = tmp;
}