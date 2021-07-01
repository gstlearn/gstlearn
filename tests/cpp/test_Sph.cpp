/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_d.h"
#include "geoslib_f.h"


#define VERBOSE 0

static double deg(double a)
{
  return(a * 180. / GV_PI);
}

static void st_test_1(void)
{
  double long1,lat1,long2,lat2,angle;
  int test;

  mestitle(1,"Testing the distance between two points");
  
  test = 1;
  while(test)
  {
    message("long1 = ");
    if (scanf("%lf",&long1) == EOF) return;
    message("lat1  = ");
    if (scanf("%lf",&lat1)  == EOF) return;
    message("long2 = ");
    if (scanf("%lf",&long2) == EOF) return;
    message("lat2  = "); 
    if (scanf("%lf",&lat2)  == EOF) return;
    angle = ut_geodetic_angular_distance(long1,lat1,
                                         long2,lat2);
    message("Angle = %lf\n",deg(angle));

    message("Continue (1) or Stop(0) : ");
    if (scanf("%d",&test) == EOF) return;
  }
}

static void st_test_2(void)
{
  double long1,lat1,long2,lat2,long3,lat3,a,b,c,A,B,C,perimeter,surface;
  double ra,rb,rc;
  int    test;

  mestitle(1,"Testing angles of a spherical triangle");

  test = 1;
  while(test)
  {
    message("long1 = ");
    if (scanf("%lf",&long1) == EOF) return;
    message("lat1  = ");
    if (scanf("%lf",&lat1)  == EOF) return;
    message("long2 = ");
    if (scanf("%lf",&long2) == EOF) return;
    message("lat2  = "); 
    if (scanf("%lf",&lat2)  == EOF) return;
    message("long3 = ");
    if (scanf("%lf",&long3) == EOF) return;
    message("lat3  = "); 
    if (scanf("%lf",&lat3)  == EOF) return;
    ut_geodetic_angles(long1,lat1,long2,lat2,long3,lat3,
                       &a,&b,&c,&A,&B,&C);
    message("a=%lf b=%lf c=%lf A=%lf B=%lf C=%lf\n",
           deg(a),deg(b),deg(c),deg(A),deg(B),deg(C));
    ra = (sin(a) == 0.) ? 1. : sin(A) / sin(a);
    rb = (sin(b) == 0.) ? 1. : sin(B) / sin(b);
    rc = (sin(c) == 0.) ? 1. : sin(C) / sin(c);
    message("ratio=%lf ratiob=%lf ratioc=%lf\n",ra,rb,rc);
    perimeter = ut_geodetic_triangle_perimeter(long1,lat1,
                                               long2,lat2,
                                               long3,lat3);
    message("Perimeter = %lf\n",deg(perimeter));
    surface = ut_geodetic_triangle_surface(long1,lat1,
                                           long2,lat2,
                                           long3,lat3);
    message("Surface = %lf\n",surface);
    
    message("Continue (1) or Stop(0) : ");
    if (scanf("%d",&test) == EOF) return;
  }
}

static void st_test_3(void)
{
  double dx,dy,s1,s2,x1,x2,y1,y2,total;
  int    nx,ny;

  mestitle(1,"Covering half-sphere with spherical triangles");

  message("nx = ");
  if (scanf("%d",&nx) == EOF) return;
  message("ny = ");
  if (scanf("%d",&ny) == EOF) return;
  dy =  90. / (double) ny;
  dx = 360. / (double) nx;
  
  total = 0.;
  for (int iy=0; iy<ny; iy++)
  {
    y1 = dy * (double) (iy);
    y2 = dy * (double) (iy+1.);
    for (int ix=0; ix<nx; ix++)
    {
      x1 = dx * (double) (ix);
      x2 = dx * (double) (ix+1.);

      s1 = ut_geodetic_triangle_surface(x1,y1,x2,y1,x1,y2);
      s2 = ut_geodetic_triangle_surface(x2,y1,x1,y2,x2,y2);
      total += s1 + s2;
    }
  }
  total /= (2. * GV_PI);
  message("Surface totale = %lf\n",total);
}

/****************************************************************************/
/*!
** Main Program
**
*****************************************************************************/
int main(int argc, char *argv[])

{
  int flag_1 = 1;
  int flag_2 = 1;
  int flag_3 = 1;
  
  /* 1.b - Connect the Geoslib Library */

  if (setup_license("Demonstration")) goto label_end;

  /* 1.c - Setup constants */

  debug_reset();
  constant_reset();
  acknowledge_Geoslib();
  constant_define("NTCAR",8);
  constant_define("NTDEC",5);
  
  if (flag_1) st_test_1();

  if (flag_2) st_test_2();
  
  if (flag_3) st_test_3();
  
label_end:
  return(0);
}

