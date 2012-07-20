/*!
  \file fpoint3d.h
  \brief addition structure for mesh_grdecl - point3d
  \author Iskhakov Ruslan
  \date 2008-05-01
 * */

#ifndef FPOINT3D_H
#define FPOINT3D_H

#include "fpoint2d.h"
#include "vector"
#define EPSILON 1.0e-7


struct tri{//нужно выделять память, там размеры массивов ниже в комментах указаны

		long int i,j,k,z,i1,j1,k1,keyF;
		int per;
		double lwt;
		double	mu;
		double	p[4];
		double	pa[3];
		double	pb[3];
		double	pc[3];
/*		double	a1;
		double	a2;
		double	a3;
		double	a4;
		double	s1;
		double	s2;
		double	s3;
		double	s4;*/
		double *points1;//4
		double *ret;//9
		double as1[12];
		long keys1;
		//Variables for mod!---------------------------------------
		double b1[3];
		double b2[3];
		double b3[3];
		double b4[3];
		double b5[3];
		double b6[3];
		double b7[3];
		double b8[3];
		double b[24];
		double *point;//4
		double *pointkey;//4
		double *localgo;//9
		double *point1;//4
		double *point2;//4
		double rast;
		double emdina;
		double local27cubs[128];//первые три столба это координаты точки пересечения, четвёртый столбец md пересечения, 5,6,7 столбцы номер куба. 8 тип точки Бля!
		long result27cubs;
		double localpoint[3];	
};
/*
void ex(double* a,double b[],int kx, int ky, int kz, double d[])//функция поиска координат куба по заданному номеру
{
	int i;
	for(i=0;i<8;i++)
	{
		a[i*pconst]=b[i*pconst]+d[0]*kx;
		a[i*pconst+1]=b[i*pconst+1]+d[1]*ky;
		a[i*pconst+2]=b[i*pconst+2]+d[2]*kz;
	}
}*/





namespace grd_ecl
{
	

  struct fpoint3d
  {
    t_double x,y,z;

    fpoint3d()
    {
      x=0.0f,y=0.0f,z=0.0f;
    }
    fpoint3d(t_double ax, t_double ay, t_double az)
    {
      x = ax;
      y = ay;
      z= az;
    }
    fpoint3d& operator += (const fpoint3d &b)
    {
      x += b.x;
      y += b.y;
      z += b.z;
      return *this;
    }
    fpoint3d& operator *= (const t_double b)
    {
      x *= b;
      y *= b;
      z *= b;
      return *this;
    }
    fpoint3d& operator *= (const fpoint3d &b)
    {
      x *= b.x;
      y *= b.y;
      z *= b.z;
      return *this;
    }
    fpoint3d& operator /= (const t_double b)
    {
      x /= b;
      y /= b;
      z /= b;
      return *this;
    }
    bool operator < (const fpoint3d &a)
    {
      if (a.z > z)
        return true;
      return false;
    }

    inline t_double get_length();

    static fpoint3d abs(const fpoint3d &a)
    {
      return fpoint3d (fabs(a.x),fabs(a.y),fabs(a.z));
    }
    
    inline void distance_to_point (const fpoint3d &a, fpoint3d &distance)
      {
        distance.x = fabs(x - a.x); 
        distance.y = fabs(y - a.y);
        distance.z = fabs(z - a.z);
      }
  //This is our place for programming
//////////////////////////////////////
//////////////////////////////////////
//////////////////////////////////////
//////////////////////////////////////
  };



  inline fpoint3d operator + (const fpoint3d &a, const fpoint3d  &b)
  {
    return fpoint3d(a.x+b.x,a.y+b.y,a.z+b.z);
  }

  inline fpoint3d operator - (const fpoint3d &a, const fpoint3d &b)
  {
    return fpoint3d(a.x-b.x,a.y-b.y,a.z-b.z);
  }

  inline fpoint3d operator / (const fpoint3d &a, const t_double d)
  {
    return fpoint3d(a.x/d,a.y/d,a.z/d);
  }

  inline t_double fpoint3d_mixed_multiplicate(const fpoint3d &a, const fpoint3d &b, const fpoint3d &c)
  {
    return a.x*b.y*c.z+a.y*b.z*c.x+a.z*b.x*c.y-a.z*b.y*c.x-a.y*b.x*c.z-a.x*b.z*c.y;
  }

  inline fpoint3d operator ^ (const fpoint3d &u, const fpoint3d &v) //vector multiplication
  {
    return fpoint3d(u.y*v.z - u.z*v.y,
                    u.z*v.x - u.x*v.z,
                    u.x*v.y - u.y*v.x);
  }

  inline fpoint3d calc_normal (const fpoint3d &a, const fpoint3d &b, const fpoint3d &c)
  {
    fpoint3d n = (b-a)^(c-a);
    t_double len = n.get_length();
    n /= len;
    return n;
  }

  inline t_double fpoint3d::get_length() //get length of radius-vector
  {
    return sqrt(x*x+y*y+z*z);
  }

  inline t_double operator *(const fpoint3d &a, const fpoint3d &b)
  {
    return a.x*b.x+a.y*b.y+a.z*b.z;
  }

  //get cosinus between plane(a,b,c) and d - normal of other plane
  inline t_double get_cosinus(const fpoint3d &a, const fpoint3d &b, const fpoint3d &c, const fpoint3d &d)
  {
    return calc_normal(a,b,c)*d;
  }

  inline t_double get_len(const fpoint3d &a, const fpoint3d &b)
  {
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
  }

  inline t_double formula_gerona3d(const fpoint3d &p1, const fpoint3d &p2, const fpoint3d &p3)
  {
    t_double a = (p2-p1).get_length();
    t_double b = (p3-p1).get_length();
    t_double c = (p3-p2).get_length();

    t_double p = (t_double)(a+b+c)/2;
    return sqrtf(p*(p-a)*(p-b)*(p-c));
  }

  inline t_double find_area_of_side(const fpoint3d &a, const fpoint3d &b, const fpoint3d &c, const fpoint3d &d)
  //get square of side (we know the side numeration)
  {
    t_double area = formula_gerona3d(a,b,c);
    area += formula_gerona3d(c,d,b);
    return area;
  }
  //calc volume of tetraidr
  inline t_double calc_tetra_volume(const fpoint3d &a, const fpoint3d &b, const fpoint3d &c, const fpoint3d &d)
  {
    t_double volume = fpoint3d_mixed_multiplicate(d-a,d-b,d-c);
    return fabs(volume/6.0f);
  }
  

  inline fpoint3d get_projection_on_all_axis_for_one_side(const boost::array <fpoint3d, 4> &vec)
  {
    fpoint3d A;

    boost::array <fpoint2d, 4> vec2;
    //YoZ
    for (size_t i = 0; i < vec.size(); i++)
      vec2[i] = fpoint2d(vec[i].y,vec[i].z);

    A.x = formula_gerona2d(vec2[0],vec2[1],vec2[2]) + formula_gerona2d(vec2[2],vec2[3],vec2[1]);

    //XoZ
    for (size_t i = 0; i < vec.size(); i++)
      vec2[i] = fpoint2d(vec[i].x,vec[i].z);

    A.y = formula_gerona2d(vec2[0],vec2[1],vec2[2]) + formula_gerona2d(vec2[2],vec2[3],vec2[1]);

    //XoY
    for (size_t i = 0; i < vec.size(); i++)
      vec2[i] = fpoint2d(vec[i].x,vec[i].y);

    A.z = formula_gerona2d(vec2[0],vec2[1],vec2[2]) + formula_gerona2d(vec2[2],vec2[3],vec2[1]);

    return A;
  }
  /*! \brief calculate intersection point for plane and segment
      \param p1,p2,p3,p4 - nodes of block edge which defines plane (currently used only first 3 points)
      \param center_side - center of side (for making decision about intersection_point inside side or not)
      \param start_point - end_point of segment
      \param delta = end_point - start_point
      \param eps_diff - difference between almost equal numbers
      \param intersection_point - return value
      \return false if no intersection*/

  inline bool
  calc_intersection_of_plane_and_segment (const fpoint3d &p1, const fpoint3d &p2,
                                          const fpoint3d &p3, const fpoint3d &p4,
                                          const fpoint3d &start_point,
                                          const fpoint3d &delta,
                                          fpoint3d &intersection_point,
                                          double eps_diff = EPSILON)
  {
    // calculate a, b, c, d for plane equation: ax+by+cz+d=0
    // a = (y1-y2)(z1+z2)+(y2-y3)(z2+z3)+(y3-y1)(z3+z1)
    // b = (z1-z2)(x1+x2)+(z2-z3)(x2+x3)+(z3-z1)(x3+x1)
    // c = (x1-x2)(y1+y2)+(x2-x3)(y2+y3)+(x3-x1)(y3+y1)
    // d = - (a * x1 + b * y1 + c * z1)
    double a,b,c,d;
    a = (p1.y-p2.y)*(p1.z+p2.z)+(p2.y-p3.y)*(p2.z+p3.z)+(p3.y-p1.y)*(p3.z+p1.z);
    b = (p1.z-p2.z)*(p1.x+p2.x)+(p2.z-p3.z)*(p2.x+p3.x)+(p3.z-p1.z)*(p3.x+p1.x);
    c = (p1.x-p2.x)*(p1.y+p2.y)+(p2.x-p3.x)*(p2.y+p3.y)+(p3.x-p1.x)*(p3.y+p1.y);
    d = -(a * p1.x + b * p1.y + c * p1.z);

    fpoint3d side_center;
    side_center.x = (p1.x + p2.x + p3.x + p4.x) / 4.;
    side_center.y = (p1.y + p2.y + p3.y + p4.y) / 4.;
    side_center.z = (p1.z + p2.z + p3.z + p4.z) / 4.;
/*
    std::cout<<"p1     ("<<p1.x<<",\t"<<p1.y<<",\t"<<p1.z<<")\n";
    std::cout<<"p2     ("<<p2.x<<",\t"<<p2.y<<",\t"<<p2.z<<")\n";
    std::cout<<"p3     ("<<p3.x<<",\t"<<p3.y<<",\t"<<p3.z<<")\n";
    std::cout<<"p4     ("<<p4.x<<",\t"<<p4.y<<",\t"<<p4.z<<")\n";
    std::cout<<"center ("<<side_center.x<<",\t"<<side_center.y<<",\t"<<side_center.z<<")\n";
    std::cout<<"start  ("<<start_point.x<<",\t"<<start_point.y<<",\t"<<start_point.z<<")\n";
    std::cout<<"delta  ("<<delta.x<<",\t"<<delta.y<<",\t"<<delta.z<<")\n";
*/

    // intersection point = start + delta * t (from parametric definition)
    // define plane vector A=[a,b,c,d] and S=[s1, s2, s3, 1]
    // S belongs to plane, so A * S = 0
    // A * start + A * delta * t = 0    =>   t = - A * start / A * delta.
    double n1 = a * start_point.x + b * start_point.y + c * start_point.z + d;
    double n2 = a * delta.x + b * delta.y + c * delta.z;
    if (!n2)
      return false;
    double t = - n1 / n2;

  // hold t in [0,1]. length of delta equal to cut length, we only reduce it.
  // So intersection_point in case t > 1 and t < 0 will not belong to plane.
    if (t > 1)
      t = 1;
    if (t < 0)
      t = 0;
  // parametric definition of line
    intersection_point.x = start_point.x + t * delta.x;
    intersection_point.y = start_point.y + t * delta.y;
    intersection_point.z = start_point.z + t * delta.z;

    // check intersection point is between p1, p2, p3, p4 (between side[0],side[1],side[2])
    //TODO: this is not precize checking, replace it with more complicated correct check
    double len = get_len (side_center, p1);
    double len_c = get_len (intersection_point, side_center);

    if (len_c < len || (fabs (len_c - len) < eps_diff))
      return true;
    else
      return false;
  }

  typedef boost::array <fpoint3d, 8>  fpoint3d_vector;

 inline double rastmd(struct tri *ara ,std::vector<fpoint3d> &v_traj,double *m)//расстояние от начала отрезка до точки пересечения
{
	//return sqrt(pow((ara->pointkey[0]-curve[long(m[0])]),2)+pow(ara->pointkey[1]-curve[long(m[0])+1],2)+pow(ara->pointkey[2]-curve[long(m[0])+2],2));
	return sqrt(pow(ara->pointkey[0]-v_traj[long(m[0])].x,2)+pow(ara->pointkey[1]-v_traj[long(m[0])].y,2)+pow(ara->pointkey[2]-v_traj[long(m[0])].z,2));
}
inline double func(double *b,double *gogo)//функция проверки совпадения предыдущей точки и текущей
{
	if ((b[0]!=gogo[4])||(b[1]!=gogo[5])||(b[2]!=gogo[6]))
	{	
		return 0; 
	}
	return 1;
}

inline double sq(double *curve,double *point,long int key)//расстояние от начала отрезка до точки пересечения.
{
	return sqrt(pow(curve[key]-point[0],2)+pow(curve[key+1]-point[1],2)+pow(curve[key+2]-point[2],2));
}

inline	double countmd(double *curve,double key,double *point,double *md,std::vector<fpoint3d> &v_traj,std::vector<t_double> &v_md)
{
	//return md[long(key/3)]+sqrt((curve[long(key)]-point[0])*(curve[long(key)]-point[0])+(curve[long(key)+1]-point[1])*(curve[long(key)+1]-point[1])+(curve[long(key)+2]-point[2])*(curve[long(key)+2]-point[2]));
	return v_md[long(key/3)]+sqrt((v_traj[long(key)].x-point[0])*(v_traj[long(key)].x-point[0])+(v_traj[long(key)].y-point[1])*(v_traj[long(key)].y-point[1])+(v_traj[long(key)].z-point[2])*(v_traj[long(key)].z-point[2]));
}

inline double plosk(double f[], std::vector<fpoint3d> &v_traj,long int key,struct tri *ara)//функция поиска пересечения отрезка с плоскостью образованной треугольником
{
	double a1,a2,a3,a4;	
	a1=(f[4]-f[1])*(f[8]-f[2])-(f[5]-f[2])*(f[7]-f[1]);
	a2=-((f[3]-f[0])*(f[8]-f[2])-(f[5]-f[2])*(f[6]-f[0]));
	a3=(f[3]-f[0])*(f[7]-f[1])-(f[4]-f[1])*(f[6]-f[0]);
	a4=-f[0]*a1-f[1]*a2-f[2]*a3;

	if((a1*(v_traj[key].x-v_traj[key+1].x)+a2*(v_traj[key].y-v_traj[key+1].y)+a3*(v_traj[key].z-v_traj[key+1].z))!=0)
	{
		return (a4+a1*(v_traj[key].x)+a2*(v_traj[key].y)+a3*(v_traj[key].z))/(a1*(v_traj[key].x-v_traj[key+1].x)+a2*(v_traj[key].y-v_traj[key+1].y)+a3*(v_traj[key].z-v_traj[key+1].z));
	}

	/*if((a1*(w[0+key]-w[3+key])+a2*(w[1+key]-w[4+key])+a3*(w[2+key]-w[5+key]))!=0)
	{
		return (a4+a1*(w[0+key])+a2*(w[1+key])+a3*(w[2+key]))/(a1*(w[0+key]-w[3+key])+a2*(w[1+key]-w[4+key])+a3*(w[2+key]-w[5+key]));
	}*/
return 5;	
}


inline double* res(std::vector<fpoint3d> &v_traj,double b[],long int key,int flag,struct tri *ara)//функция поиска точки пересечения отрезка с треугольником
{
	 double s1,s2,s3,s4,a1,a2,a3,a4,mu;
	mu=plosk(b,v_traj,key,ara);
	if(((mu>=0)&&(mu<=1)))
		ara->p[0]=v_traj[key].x+mu*(v_traj[key+1].x-v_traj[key].x);
		ara->p[1]=v_traj[key].y+mu*(v_traj[key+1].y-v_traj[key].y);
		ara->p[2]=v_traj[key].z+mu*(v_traj[key+1].z-v_traj[key].z);
	{
	//	ara->p[0]=curve[0+key]+mu*(curve[3+key]-curve[0+key]);
	//	ara->p[1]=curve[1+key]+mu*(curve[4+key]-curve[1+key]);
	//	ara->p[2]=curve[2+key]+mu*(curve[5+key]-curve[2+key]);
		if((flag==1)||(flag==2))
		{
			a1=b[0]-b[3];
			a2=b[1]-b[4];
			a3=ara->p[0]-b[3];
			a4=ara->p[1]-b[4]; //получили вектора A и B

			s1=0.5*fabs(a1*a4-a2*a3);

			a1=b[6]-b[0];
			a2=b[7]-b[1];
			a3=ara->p[0]-b[0];
			a4=ara->p[1]-b[1]; //получили вектора A и B

			s2=0.5*fabs(a1*a4-a2*a3);

			a1=b[6]-b[3];
			a2=b[7]-b[4];
			a3=ara->p[0]-b[3];
			a4=ara->p[1]-b[4]; //получили вектора A и B

			s3=0.5*fabs(a1*a4-a2*a3);

			a1=b[0]-b[3];
			a2=b[1]-b[4];
			a3=b[6]-b[3];
			a4=b[7]-b[4]; //получили вектора A и B
			
			s4=0.5*fabs(a1*a4-a2*a3);

			if(((s1==0)||(s2==0)||(s3==0))&&(s4==s1+s2+s3))
			{
				ara->p[3]=0;//пересечение с ребром
				return ara->p;
			}

			if(s4==s1+s2+s3)
			{
				ara->p[3]=1;//пересечение с гранью
				return ara->p;
			}
		}

			if((flag==3)||(flag==4))
		{
			a1=b[1]-b[4];
			a2=b[2]-b[5];
			a3=ara->p[1]-b[4];
			a4=ara->p[2]-b[5]; //получили вектора A и B

			s1=0.5*fabs((a1)*(a4)-(a2)*(a3));

			a1=b[7]-b[1];
			a2=b[8]-b[2];
			a3=ara->p[1]-b[1];
			a4=ara->p[2]-b[2]; //получили вектора A и B

			s2=0.5*fabs(a1*a4-a2*a3);

			a1=b[7]-b[4];
			a2=b[8]-b[5];
			a3=ara->p[1]-b[4];
			a4=ara->p[2]-b[5]; //получили вектора A и B

			s3=0.5*fabs(a1*a4-a2*a3);

			a1=b[1]-b[4];
			a2=b[2]-b[5];
			a3=b[7]-b[4];
			a4=b[8]-b[5]; //получили вектора A и B
			
			s4=0.5*fabs(a1*a4-a2*a3);
			if(((s1==0)||(s2==0)||(s3==0))&&(s4==s1+s2+s3))
			{
				ara->p[3]=0;//пересечение с ребром
				return ara->p;
			}

			if(s4==s1+s2+s3)
			{
				ara->p[3]=1;
				return ara->p;
			}
		}
		

		if((flag==5)||(flag==6))
		{
			a1=b[0]-b[3];
			a2=b[2]-b[5];
			a3=ara->p[0]-b[3];
			a4=ara->p[2]-b[5]; //получили вектора A и B

			s1=0.5*fabs(a1*a4-a2*a3);

			a1=b[6]-b[0];
			a2=b[8]-b[2];
			a3=ara->p[0]-b[0];
			a4=ara->p[2]-b[2]; //получили вектора A и B

			s2=0.5*fabs(a1*a4-a2*a3);

			a1=b[6]-b[3];
			a2=b[8]-b[5];
			a3=ara->p[0]-b[3];
			a4=ara->p[2]-b[5]; //получили вектора A и B

			s3=0.5*fabs(a1*a4-a2*a3);

			a1=b[0]-b[3];
			a2=b[2]-b[5];
			a3=b[6]-b[3];
			a4=b[8]-b[5]; //получили вектора A и B
			
			s4=0.5*fabs(a1*a4-a2*a3);
			if(((s1==0)||(s2==0)||(s3==0))&&(s4==s1+s2+s3))
			{
				ara->p[3]=0;//пересечение с ребром
				return ara->p;
			}

			if(s4==s1+s2+s3)
			{
				ara->p[3]=1;
				return ara->p;
			}
		}
	}
	return NULL;
}





//Основной алгоритм

inline void cpypoint(double *points,double *point,long int &flag,double *gogo)//просто функция копирования
{
	points[flag]=point[0];
	points[flag+1]=point[1];
	points[flag+2]=point[2];
//points[flag]=gogo[1];
//points[flag+1]=gogo[2];
//points[flag+2]=gogo[3];
}
inline void cpygo(double *points,double *gogo,long int &flag)//просто функция копирования
{
	points[flag]=gogo[4];
	points[flag+1]=gogo[5];
	points[flag+2]=gogo[6];
}


inline long inf27points(double massiv[],long key, double md)//функция ищет из найденных точек пересечения с 27-ю кубами самую ближнюю
{
	double B;
	long k,i;
	B=massiv[3];
	k=0;
	for(i=3;i<=key*8-4;i+=8)
	{
		if ((B >= massiv[i])&&(B>md))
		{
			k=i;
			B=massiv[i];
		}
	}

	k=int(k/8);
	return k;
}


inline double func1(double *points,long int flag, struct tri *ara,std::vector<fpoint3d> v_traj, double *gogo)//функция проверяющая нашли мы правильную точку или нет.
{

				//ara->emdina=sqrt(pow((ara->point1[0]-curve[long(gogo[0])]),2)+pow((ara->point1[1]-curve[long(gogo[0])+1]),2)+pow((ara->point1[2]-curve[long(gogo[0])+2]),2));
				ara->emdina=sqrt(pow((ara->point1[0]-v_traj[long(gogo[0])].x),2)+pow((ara->point1[1]-v_traj[long(gogo[0])].y),2)+pow((ara->point1[2]-v_traj[long(gogo[0])].z),2));
				//if((ara->localpoint[0]==ara->point1[0])&&(ara->localpoint[1]==ara->point1[1])&&(ara->localpoint[2]==ara->point1[2]))
				//return 0;
				if(ara->rast!=0)
				{
				if((ara->emdina)>=(ara->rast))
				{
					return 1;//если правильная
				}
			return 0;//если нет
					}
					else if(ara->emdina>ara->localpoint[0])
							{
								return 1;
							}
				return 0;
}
};//namespace grd_ecl
#endif
