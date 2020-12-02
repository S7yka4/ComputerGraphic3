#define _USE_MATH_DEFINES
#include <math.h>

#include "Render.h"
#include "Math.h"
#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>

double t_max = 0;

void DrawSquare(double v1[], double v2[], double v3[], double v4[], int flag = 0, int i = 0, int n = 0)
{

	glBegin(GL_QUADS);
	glVertex3dv(v1);
	glVertex3dv(v2);
	glVertex3dv(v3);
	glVertex3dv(v4);
	glEnd();
}
void DrawWall(double v1[], double v2[], int flag = 0, int i = 0, int n = 0)
{
	double tmp1[3];
	double tmp2[3];
	tmp1[0] = v1[0];
	tmp1[1] = v1[1];
	tmp1[2] = v1[2] - 1;
	tmp2[0] = v2[0];
	tmp2[1] = v2[1];
	tmp2[2] = v2[2] - 1;
	//glNormal3d(CalculateNorm(v1, tmp1, v2)[0], CalculateNorm(v1, tmp1, v2)[1], CalculateNorm(v1, tmp1, v2)[2]);
	DrawSquare(v1, v2, tmp2, tmp1, flag, i, n);
}
void DrawCube(double P1[], double P2[], double P3[], double P4[])
{
	double PA[] = { P1[0] + 0.5,P1[1] + 0.5,P1[2] - 0.5 };
	double PB[] = { P1[0] + 0.5,P1[1] - 0.5,P1[2] - 0.5 };
	double PC[] = { P1[0] - 0.5,P1[1] - 0.5,P1[2] - 0.5 };
	double PD[] = { P1[0] - 0.5,P1[1] + 0.5,P1[2] - 0.5 };
	/*double PA[] = { 0.5, 0.5, -0.5 };
	double PB[] = { 0.5, -0.5,-0.5 };
	double PC[] = { -0.5, -0.5, -0.5 };
	double PD[] = { -0.5, 0.5, -0.5 };*/
	DrawSquare(PA, PB, PC, PD);
	double PA1[] = { P1[0] + 0.5,P1[1] + 0.5,P1[2] + 0.5 };
	double PB1[] = { P1[0] + 0.5,P1[1] - 0.5,P1[2] + 0.5 };
	double PC1[] = { P1[0] - 0.5,P1[1] - 0.5,P1[2] + 0.5 };
	double PD1[] = { P1[0] - 0.5,P1[1] + 0.5,P1[2] + 0.5 };
	DrawSquare(PA1, PB1, PC1, PD1);
	DrawWall(PA1, PB1);
	DrawWall(PC1, PB1);
	DrawWall(PC1, PD1);
	DrawWall(PD1, PA1);
}

 double Bez(double p0,double p1, double p2, double p3, double t)
{
	
	 return ((1 - t) * (1 - t) * (1 - t) * p0) + (3 * t*(1 - t) * (1 - t) * p1) + (3 * t * t * (1 - t)*p2) + (t * t * t * p3);
} 
 void  DrawBaseLines(double PointA[], double PointB[], double PointC[], double PointD[])
 {
	 glColor3d(1, 0, 0);
	 glBegin(GL_LINE_STRIP);
	 glVertex3dv(PointA);
	 glVertex3dv(PointB);

	 glVertex3dv(PointD);
	 glVertex3dv(PointC);
	 glEnd();
	 

 }
 void DrawBez(double PointA[], double PointB[], double PointC[], double PointD[])
 {
	 double P[3];

	 DrawBaseLines(PointA, PointB, PointC, PointD);
	 glColor3d(1, 0, 0);
	 glBegin(GL_LINE_STRIP);
	 for (double t = 0; t < 1.0001; t = t + 0.01)
	 {
		 P[0] = Bez(PointA[0], PointB[0], PointC[0], PointD[0], t);
		 P[1] = Bez(PointA[1], PointB[1], PointC[1], PointD[1], t);
		 P[2] = Bez(PointA[2], PointB[2], PointC[2], PointD[2], t);
		 glVertex3dv(P);
	 }
	 glEnd();
 }

 void MoveCube(double PointA[], double PointB[], double PointC[], double PointD[], double delta_time)
 {
	 double P[3];
	 double oldP[] = { PointA[0],PointA[1],PointA[2] };
	 static bool flagReverse = false;

	 if (!flagReverse)
	 {
		 t_max += delta_time / 20; //t_max становится = 1 за 5 секунд
		 if (t_max > 1)
		 {
			 t_max = 1; //после обнуляется
			 flagReverse = !flagReverse;
		 }
	 }
	 else
	 {
		 t_max -= delta_time / 20; //t_max становится = 1 за 5 секунд
		 if (t_max < 0)
		 {
			 t_max = 0; //после обнуляется
			 flagReverse = !flagReverse;
		 }
	 }
	 
			 glPushMatrix();
			 for (double t = 0; t <= t_max; t = t + 0.001)
			 {
				 P[0] = Bez(PointA[0], PointB[0], PointC[0], PointD[0], t);
				 P[1] = Bez(PointA[1], PointB[1], PointC[1], PointD[1], t);
				 P[2] = Bez(PointA[2], PointB[2], PointC[2], PointD[2], t);
				 //glRotated(acos((P[2] - oldP[2]) / ((sqrt(pow(P[0] - oldP[0], 2) + pow(P[1] - oldP[1], 2) + pow(P[2] - oldP[2], 2)))))*180/M_PI,P[0] - oldP[0], P[1] - oldP[1], P[2] - oldP[2]);
				 glTranslated(P[0] - oldP[0], P[1] - oldP[1], P[2] - oldP[2]);
				 
				 oldP[0] = P[0];
				 oldP[1] = P[1];
				 oldP[2] = P[2];
			 }
			 DrawCube(PointA,PointB,PointC, PointD);
			 glPopMatrix();
 }



 double* Erm(double VecA[], double VecA1[]  , double VecC[] , double VecC1[], double t)
 {
	 double Result[3];
	 for(int i=0;i<3;i++)

         Result[i]=VecA[i] * (2 * t * t * t - 3 * t * t + 1) + VecC1[i]*(-2 * t * t * t + 3 * t * t) + VecA1[i] * (t * t * t - 2 * t * t + t) + VecC[i]*(t * t * t - t * t);
	 return Result;
 }
 void DrawBaseVecs(double PointA[], double PointB[], double PointC[], double PointD[])
 {
	 glColor3d(0, 1, 0);
	 glBegin(GL_LINES);
	 glVertex3dv(PointA);
	 glVertex3dv(PointB);
	 glVertex3dv(PointC);
	 glVertex3dv(PointD);
	 glEnd();
	 
	 glColor3d(1, 0, 1);
	 glBegin(GL_POINTS);
	 glVertex3dv(PointA);
	 glVertex3dv(PointC);
	 glEnd();
 }
 void DrawErm(double PointA[], double PointB[], double PointC[], double PointD[])
 {
	 double P[3];
	 double VecAB[] = {PointB[0]-PointA[0], PointB[1]-PointA[1],PointB[2] - PointA[2] };
	 double VecCD[] = { PointD[0] - PointC[0], PointD[1] - PointC[1],PointD[2] - PointC[2] };
	 DrawBaseVecs(PointA, PointB, PointC, PointD);
	 glColor3d(0, 1, 0);
	 glBegin(GL_LINE_STRIP);
	 for (double i = 0; i < 1.0001; i = i + 0.01)
	 {
		 P[0] = *(Erm(PointA, VecAB, VecCD,   PointC, i));
		 P[1] = *(Erm(PointA, VecAB, VecCD,  PointC, i) +1);
		 P[2] = *(Erm(PointA, VecAB, VecCD,  PointC, i) +2);
		 
		 glVertex3dv(P);
	 }
	 glEnd();


	
 }


 /*Зеленое - все что связанно с кривой эрмита, красное все что связано с кривой безье*/

void Render(double delta_time)
{

	double PointA1[] = { 2,1.5, 0  };
	double PointB1[] = { 7, 3.5,0 };
	double PointC1[] = { 1.5, -0.5, 0 };
	double PointD1[] = { 6, -2.5,0 };

	
	DrawErm(PointA1, PointB1, PointC1, PointD1);//рисует кривую эрмиат и вектора
	DrawBez(PointA1, PointB1, PointC1, PointD1);//рисует кривую безье и отрезки
		MoveCube(PointA1, PointB1, PointC1, PointD1, delta_time);//двигает куб по кривой безье 


};



