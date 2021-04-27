#define _CRT_SECURE_NO_WARNINGS
#define MAXIMUM_DEG 165
#define MINIMUM_DEG	110
#define SHOOTING_FORCE	1.1f
#define EXPLOSION_TIME	1
#include "glut.h"
#include <math.h>
#include <time.h>
#include <vector>
#include <windows.h>

using namespace std;


const int W = 600;
const int H = 600;
const double PI = 3.14;
const int GSZ = 150;

// texture matrices
const int  TSZ = 512;
unsigned char tx0[TSZ][TSZ][3]; // 3 stands for RGB
unsigned char tx1[1024][1024][3]; // 3 stands for RGB
unsigned char tx2[512][1024][3]; // 3 stands for RGB
unsigned char tx3[512][1024][3]; // tank
unsigned char tx4[512][1024][3]; // toptree 1
unsigned char tx5[512][1024][3]; // bottomtree 1
unsigned char tx6[512][1024][3]; // totah
unsigned char tx7[512][1024][3]; // explotion texture
unsigned char tx8[512][1024][3]; // chain back&front
unsigned char tx9[512][1024][3]; // chain rest

//trees
struct point
{
	int x;
	int z;
};
struct point Trees[100];
/// 

double ground[GSZ][GSZ] = { 0 };
double red[GSZ][GSZ], green[GSZ][GSZ], blue[GSZ][GSZ];
double offset = 0;

unsigned char* bmp; // array of pixels

double eyex = -1, eyey = 19, eyez = 15;

double dx = 0, dy = 0, dz = 0;
// ego-motion
double speed = 0, angular_speed = 0;
double sight_angle = PI; // initial sight
double pitch = -0.2;
double dir[3]; // x,y,z
// airplane
double aspeed = 0, aangular_speed = 0;
double yaw = PI;
double apitch1v = 0;
double apitch2h = 0;
double adir[3] = { sin(yaw),sin(apitch1v),cos(yaw) };
double ax = 0, ay = 12, az = 0;

bool shouldShoot = false;
bool shouldExplode = false;
bool is_building = true;
double t=0 , tExplostion =0;
double a = -20.8;
double yawi;
double  verticalApitch;
double horazintalApitch;
double xi, yi, zi;
double vx=50, vy = 15, vz = 50;
double xf=0, yf=0, zf=0;

const double l1 = 0.064;
const double l2 = 0.11;

double shootTime = 0;
double explosionStartTime = 0;
void UpdateTerrain1();
void UpdateTerrain2();
void Smooth();

void ReadPicture(char* fname)
{
	FILE* pf;
	BITMAPFILEHEADER bf;
	BITMAPINFOHEADER bi;
	pf = fopen(fname, "rb");

	fread(&bf, sizeof(bf), 1, pf);
	fread(&bi, sizeof(bi), 1, pf);
	bmp = (unsigned char*)malloc(3 * bi.biWidth * bi.biHeight); // 3 stands for BGR
	fread(bmp, 1, 3 * bi.biWidth * bi.biHeight, pf);
	fclose(pf);
}

void SetTexture(int txnum)
{
	int i, j, k;
	switch (txnum)
	{
	case 0: // brick wall
		for (i = 0; i < TSZ; i++)
			for (j = 0; j < TSZ; j++)
			{
				k = rand() % 20;
				if (i % (TSZ / 2) < 20 ||  // horizontal separation
					i > TSZ / 2 && j % (TSZ / 2) < 20 || // upper part of bricks
					i <= TSZ / 2 && (j > TSZ / 4 && j < TSZ / 4 + 20 || j>3 * TSZ / 4 && j < 3 * TSZ / 4 + 20)  // lower paert of bricks
					)
				{    // separating layer
					tx0[i][j][0] = 150 + k; //red
					tx0[i][j][1] = 150 + k; //green 
					tx0[i][j][2] = 150 + k; //blue
				}
				else // brick
				{
					tx0[i][j][0] = 150 + k; //red
					tx0[i][j][1] = 50 + k; //green 
					tx0[i][j][2] = k; //blue
				}
			}
		break;
	case 1:
		for (i = 0; i < 1024; i++)
			for (j = 0; j < 1024; j++)
			{
				tx1[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx1[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx1[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	case 2:
		for (i = 0; i < 512; i++)
			for (j = 0; j < 1024; j++)
			{
				tx2[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx2[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx2[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;

	case 3:
		for (i = 0; i < 512; i++)
			for (j = 0; j < 1024; j++)
			{
				tx3[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx3[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx3[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	case 4:
		for (i = 0; i < 512; i++)
			for (j = 0; j < 1024; j++)
			{
				tx4[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx4[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx4[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	case 5:
		for (i = 0; i < 512; i++)
			for (j = 0; j < 1024; j++)
			{
				tx5[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx5[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx5[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	case 6:
		for (i = 0; i < 512; i++)
			for (j = 0; j < 1024; j++)
			{
				tx6[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx6[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx6[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	case 7:
		for (i = 0; i < 512; i++)
			for (j = 0; j < 1024; j++)
			{
				tx7[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx7[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx7[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	case 8:
		for (i = 0; i < 512; i++)
			for (j = 0; j < 1024; j++)
			{
				tx8[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx8[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx8[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	case 9:
		for (i = 0; i < 512; i++)
			for (j = 0; j < 1024; j++)
			{
				tx9[i][j][0] = bmp[(i * 1024 + j) * 3 + 2];  //red
				tx9[i][j][1] = bmp[(i * 1024 + j) * 3 + 1];// green
				tx9[i][j][2] = bmp[(i * 1024 + j) * 3];//blue
			}
		break;
	}
}
// topr and bottom r are the top and bottom radiuses
void DrawTexCylinder1(int n, double topr, double bottomr, int tnum, double hrepeats)
{
	double alpha, teta = 2 * PI / n;
	int side_counter;
	double tex_part = hrepeats / n; // part of a texture that covers one side
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, tnum); // texture #0
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);// GL_REPLACE);// GL_MODULATE); 

	for (alpha = 0, side_counter = 0; alpha <= 2 * PI; alpha += teta, side_counter++)
	{
		glBegin(GL_POLYGON);
		glColor3d(0.5 + fabs(sin(alpha)) / 2, 0.5 + fabs(sin(alpha)) / 2, 0.5 + fabs(sin(alpha) / 2));
		glTexCoord2d(side_counter * tex_part, 1);	glVertex3d(topr * sin(alpha), 1, topr * cos(alpha)); //1
//		glColor3d(0.2 + fabs(sin(alpha + teta) / 2), 0.2 + fabs(sin(alpha + teta) / 2), 0.2 +fabs( sin(alpha + teta) / 2));
		glTexCoord2d((side_counter + 1) * tex_part, 1);		glVertex3d(topr * sin(alpha + teta), 1, topr * cos(alpha + teta));//2
		glTexCoord2d((side_counter + 1) * tex_part, 0);		glVertex3d(bottomr * sin(alpha + teta), 0, bottomr * cos(alpha + teta));//3
//		glColor3d(0.2 + fabs(sin(alpha) / 2), 0.2 +fabs( sin(alpha) / 2), 0.2 +fabs( sin(alpha) / 2));
		glTexCoord2d(side_counter * tex_part, 0);		glVertex3d(bottomr * sin(alpha), 0, bottomr * cos(alpha));// 4
		glEnd();
	}

	glDisable(GL_TEXTURE_2D);

}
// topr and bottom r are the top and bottom radiuses
void DrawTexCylinder2(int n, double topr, double bottomr, int tnum, double hrepeats,
	double txbottom, double txtop)
{
	double alpha, teta = 2 * PI / n;
	int side_counter;
	double tex_part = hrepeats / n; // part of a texture that covers one side
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, tnum); // texture #0
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);// GL_MODULATE); 

	for (alpha = 0, side_counter = 0; alpha <= 2 * PI; alpha += teta, side_counter++)
	{
		glBegin(GL_POLYGON);
		glColor3d(0.5 + fabs(sin(alpha)) / 2, 0.5 + fabs(sin(alpha)) / 2, 0.5 + fabs(sin(alpha) / 2));
		glTexCoord2d(side_counter * tex_part, txtop);	glVertex3d(topr * sin(alpha), 1, topr * cos(alpha)); //1
//		glColor3d(0.2 + fabs(sin(alpha + teta) / 2), 0.2 + fabs(sin(alpha + teta) / 2), 0.2 +fabs( sin(alpha + teta) / 2));
		glTexCoord2d((side_counter + 1) * tex_part, txtop);		glVertex3d(topr * sin(alpha + teta), 1, topr * cos(alpha + teta));//2
		glTexCoord2d((side_counter + 1) * tex_part, txbottom);		glVertex3d(bottomr * sin(alpha + teta), 0, bottomr * cos(alpha + teta));//3
//		glColor3d(0.2 + fabs(sin(alpha) / 2), 0.2 +fabs( sin(alpha) / 2), 0.2 +fabs( sin(alpha) / 2));
		glTexCoord2d(side_counter * tex_part, txbottom);		glVertex3d(bottomr * sin(alpha), 0, bottomr * cos(alpha));// 4
		glEnd();
	}
	glDisable(GL_TEXTURE_2D);
}

void DrawTexSphere(int slices, int stacks, int tnum, double hrep, double vrep)
{

	double beta = PI / slices, gamma; // gamma is running angle from -90 to 90
	int i;
	double vpart = vrep / stacks;

	for (gamma = -PI / 2, i = 0; gamma <= PI / 2; gamma += beta, i++)
	{
		glPushMatrix();
		//		glRotated(i * offset/10, 0, 1, 0);
		glTranslated(0, sin(gamma), 0);
		glScaled(1, sin(gamma + beta) - sin(gamma), 1);
		DrawTexCylinder2(stacks, cos(gamma + beta), cos(gamma), tnum, hrep, vpart * i, vpart * (i + 1));
		glPopMatrix();
	}
}

void init()
{
	int i, j;
	srand(time(0));
	//red  green  blue
	glClearColor(0.0, 0.3, 0.6, 0);// color of window background
	glEnable(GL_DEPTH_TEST); // show the nearest object

	for (i = 0; i < GSZ; i++)
		for (j = 0; j < GSZ; j++)
			ground[i][j] = 0.1;//3*sin(j/3.0)+2*sin(i/2.0);
	for (i = 1; i <= 1000; i++)
		UpdateTerrain1();
	for (i = 1; i <= 1200; i++)
		UpdateTerrain2();
	Smooth();
	Smooth();
	
	for (int i = 0; i <= 40; i++) {
	int x = rand() % (GSZ - 2);
	int z = rand() % (GSZ - 2);

		while (ground[z][x] <= 1) {
			x = rand() % (GSZ - 2);
			z = rand() % (GSZ - 2);
	
		}
		Trees[i].x = x - GSZ / 2;
		Trees[i].z = z - GSZ / 2;
	}
	

	glEnable(GL_NORMALIZE); // used for lighting

	char name[30];
	strcpy(name, "sky.bmp");
	ReadPicture(name);
	SetTexture(2); // sky
	glBindTexture(GL_TEXTURE_2D, 2); // setting texture #2
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, tx2);

	strcpy(name, "tankpic.bmp");
	ReadPicture(name);
	SetTexture(3); // sky
	glBindTexture(GL_TEXTURE_2D, 3); // setting texture #2
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, tx3);

	strcpy(name, "toptree1.bmp");
	ReadPicture(name);
	SetTexture(4); // top tree1
	glBindTexture(GL_TEXTURE_2D, 4); // setting texture #4
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, tx4);

	strcpy(name, "bottomtree1.bmp");
	ReadPicture(name);
	SetTexture(5); // bottom tree1
	glBindTexture(GL_TEXTURE_2D, 5); // setting texture #4
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, tx5);

	strcpy(name, "totah.bmp");
	ReadPicture(name);
	SetTexture(6);
	glBindTexture(GL_TEXTURE_2D, 6); // setting texture #4
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, tx6);

	strcpy(name, "exp.bmp");
	ReadPicture(name);
	SetTexture(7);
	glBindTexture(GL_TEXTURE_2D, 7); // setting texture #4
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, tx7);

	strcpy(name, "chain.bmp");
	ReadPicture(name);
	SetTexture(8);
	glBindTexture(GL_TEXTURE_2D, 8); // setting texture #4
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, tx8);

	strcpy(name, "chain2.bmp");
	ReadPicture(name);
	SetTexture(9);
	glBindTexture(GL_TEXTURE_2D, 9); // setting texture #4
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // horizontal 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // vertical
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1024, 512, 0, GL_RGB, GL_UNSIGNED_BYTE, tx9);
}

void UpdateTerrain()
{
	int i, j;
	double delta = 0.05;

	for (i = 0; i < GSZ; i++)
		for (j = 0; j < GSZ; j++)
		{
			if (rand() % 100 >= 50)
				delta = -delta;
			ground[i][j] += delta;
		}
}

void UpdateTerrain1()
{
	int i, j;
	double delta = 0.05, a, b;
	int x1, y1, x2, y2;

	if (rand() % 2 == 0)
		delta = -delta;

	x1 = rand() % GSZ;
	y1 = rand() % GSZ;
	x2 = rand() % GSZ;
	y2 = rand() % GSZ;

	if (x1 == x2) // exception
		return;

	a = (y2 - y1) / (double)(x2 - x1);
	b = y1 - a * x1;

	for (i = 0; i < GSZ; i++)
		for (j = 0; j < GSZ; j++)
		{
			if (i < j * a + b)
				ground[i][j] += delta;
			else
				ground[i][j] -= delta;
		}
}
// random walk
void UpdateTerrain2()
{
	int steps = 5000;
	double delta = 0.01;
	int x, y;
	if (rand() % 2 == 0)
		delta = -delta;

	x = rand() % GSZ;
	y = rand() % GSZ;

	for (int counter = 1; counter <= steps; counter++)
	{
		ground[y][x] += delta;
		// choose new direction
		switch (rand() % 4)
		{
		case 0: // up
			y++;
			break;
		case 1: // down
			y--;
			break;
		case 2: // right
			x++;
			break;
		case 3: // left
			x--;
			break;
		}
		// in case that we go out of matrix we use the formula:
		x = (x + GSZ) % GSZ;
		y = (y + GSZ) % GSZ;
	}
}
void Smooth()
{
	int i, j;
	double tmp[GSZ][GSZ];

	for (i = 1; i < GSZ - 1; i++)
		for (j = 1; j < GSZ - 1; j++)
			tmp[i][j] = (ground[i - 1][j - 1] + 2 * ground[i - 1][j] + ground[i - 1][j + 1] +
				2 * ground[i][j - 1] + 4 * ground[i][j] + 2 * ground[i][j + 1] +
				ground[i + 1][j - 1] + 2 * ground[i + 1][j] + ground[i + 1][j + 1]) / 16.0;
	// copy tmp to ground
	for (i = 1; i < GSZ - 1; i++)
		for (j = 1; j < GSZ - 1; j++)
			ground[i][j] = tmp[i][j];
}
void SetColor(double h)
{

	//	glColor3d(fabs(h / 3), (h+3)/6, fabs(sin(h)));
	h = fabs(h);
	if (h < 0.15) // sand
		glColor3d(0.8, 0.7, 0.5);
	else
		if (h < 3) // grass
			glColor3d(0.1 + h / 11, 0.5 - h / 10, 0);
		else // snow
			glColor3d(h / 5, h / 5, h / 4);

}
void SetNormal(int row, int col)
{
	// if  it was defined glEnable(GL_NORMALIZE) then we don't need
	// to normalize it manually
	glNormal3d(ground[row][col + 1] - ground[row][col], 1,
		ground[row + 1][col] - ground[row][col]);
}
void DrawLitFloor()
{
	int i, j;

	glColor3d(0.8, 0.8, 1); // blue

	for (i = 1; i < GSZ - 2; i++)
		for (j = 1; j < GSZ - 2; j++)
		{
			glBegin(GL_POLYGON); // lines parallel to X
			SetNormal(i, j);
			glVertex3d(j - GSZ / 2, ground[i][j], i - GSZ / 2);
			SetNormal(i - 1, j);
			glVertex3d(j - GSZ / 2, ground[i - 1][j], i - 1 - GSZ / 2);
			SetNormal(i - 1, j - 1);
			glVertex3d(j - 1 - GSZ / 2, ground[i - 1][j - 1], i - 1 - GSZ / 2);
			SetNormal(i, j - 1);
			glVertex3d(j - 1 - GSZ / 2, ground[i][j - 1], i - GSZ / 2);
			glEnd();
		}


	// water
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4d(0.2, 0.4, 0.6, 0.7);
	glBegin(GL_POLYGON);
	glVertex3d(-GSZ / 2, 0, -GSZ / 2);
	glVertex3d(-GSZ / 2, 0, GSZ / 2);
	glVertex3d(GSZ / 2, 0, GSZ / 2);
	glVertex3d(GSZ / 2, 0, -GSZ / 2);
	glEnd();
	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);


}
void DrawFloor()
{
	int i, j;
	glColor3d(0.8, 0.8, 0.8); // blue
	for (i = 1; i < GSZ - 2; i++)
		for (j = 1; j < GSZ - 2; j++)
		{
			glBegin(GL_POLYGON); // lines parallel to X
			SetColor(ground[i][j]);
			//			SetNormal(i, j);
			glVertex3d(j - GSZ / 2, ground[i][j], i - GSZ / 2);
			SetColor(ground[i - 1][j]);
			//			SetNormal(i - 1, j);
			glVertex3d(j - GSZ / 2, ground[i - 1][j], i - 1 - GSZ / 2);
			SetColor(ground[i - 1][j - 1]);
			//			SetNormal(i - 1, j - 1);
			glVertex3d(j - 1 - GSZ / 2, ground[i - 1][j - 1], i - 1 - GSZ / 2);
			SetColor(ground[i][j - 1]);
			//			SetNormal(i, j - 1);
			glVertex3d(j - 1 - GSZ / 2, ground[i][j - 1], i - GSZ / 2);
			glEnd();
		}

	// water
//	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4d(0.2, 0.4, 0.6, 0.7);
	glBegin(GL_POLYGON);
	glVertex3d(-GSZ / 2, 0, -GSZ / 2);
	glVertex3d(-GSZ / 2, 0, GSZ / 2);
	glVertex3d(GSZ / 2, 0, GSZ / 2);
	glVertex3d(GSZ / 2, 0, -GSZ / 2);
	glEnd();
	glDisable(GL_BLEND);
	//	glEnable(GL_LIGHTING);
}
void DrawAxes()
{
	glLineWidth(2);

	// x
	glColor3d(1, 0, 0);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(12, 0, 0);
	glEnd();
	// y
	glColor3d(0, 1, 0);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 12, 0);
	glEnd();
	// z
	glColor3d(1, 1, 0);
	glBegin(GL_LINES);
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 12);
	glEnd();

	glLineWidth(1);

}
void DrawColorCube()
{
	// top face
	glBegin(GL_POLYGON);
	glColor3d(1, 0, 0); // red
	glVertex3d(1, 1, 1);
	glColor3d(1, 1, 0); // yellow
	glVertex3d(1, 1, -1);
	glColor3d(0, 1, 0); // green
	glVertex3d(-1, 1, -1);
	glColor3d(1, 1, 1); // white
	glVertex3d(-1, 1, 1);
	glEnd();
	// bottom face
	glBegin(GL_POLYGON);
	glColor3d(1, 0, 1); // magenta
	glVertex3d(1, -1, 1);
	glColor3d(0, 0, 0); // black
	glVertex3d(1, -1, -1);
	glColor3d(0, 1, 1); // cyan
	glVertex3d(-1, -1, -1);
	glColor3d(0, 0, 1); // blue
	glVertex3d(-1, -1, 1);
	glEnd();
	// left face
	glBegin(GL_POLYGON);
	glColor3d(1, 1, 1); // white
	glVertex3d(-1, 1, 1);
	glColor3d(0, 1, 0); // green
	glVertex3d(-1, 1, -1);
	glColor3d(0, 1, 1); // cyan
	glVertex3d(-1, -1, -1);
	glColor3d(0, 0, 1); // blue
	glVertex3d(-1, -1, 1);
	glEnd();
	// right face
	glBegin(GL_POLYGON);
	glColor3d(1, 0, 0); // red
	glVertex3d(1, 1, 1);
	glColor3d(1, 1, 0); // yellow
	glVertex3d(1, 1, -1);
	glColor3d(0, 0, 0); // black
	glVertex3d(1, -1, -1);
	glColor3d(1, 0, 1); // magenta
	glVertex3d(1, -1, 1);
	glEnd();
	// back
	glBegin(GL_POLYGON);
	glColor3d(1, 1, 0); // yellow
	glVertex3d(1, 1, -1);
	glColor3d(0, 1, 0); // green
	glVertex3d(-1, 1, -1);
	glColor3d(0, 1, 1); // cyan
	glVertex3d(-1, -1, -1);
	glColor3d(0, 0, 0); // black
	glVertex3d(1, -1, -1);
	glEnd();

	// front face
	glBegin(GL_POLYGON);
	glColor3d(1, 1, 1); // white
	glVertex3d(-1, 1, 1);
	glColor3d(1, 0, 0); // red
	glVertex3d(1, 1, 1);
	glColor3d(1, 0, 1); // magenta
	glVertex3d(1, -1, 1);
	glColor3d(0, 0, 1); // blue
	glVertex3d(-1, -1, 1);
	glEnd();


}
// topr and bottom r are the top and bottom radiuses
void DrawCylinder1(int n, double topr, double bottomr)
{
	double alpha, teta = 2 * PI / n;
	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		glBegin(GL_POLYGON);
		glColor3d(fabs(sin(alpha) / 2), fabs(sin(alpha)), fabs(cos(alpha)) / 2);
		glVertex3d(topr * sin(alpha), 1, topr * cos(alpha)); //1

		glVertex3d(topr * sin(alpha + teta), 1, topr * cos(alpha + teta));//2
		glColor3d(fabs(sin(alpha) / 2), fabs(sin(alpha)), fabs(cos(alpha)) / 2);
		glVertex3d(bottomr * sin(alpha + teta), 0, bottomr * cos(alpha + teta));//3
		glVertex3d(bottomr * sin(alpha), 0, bottomr * cos(alpha));// 4
		glEnd();
	}
}
// topr and bottom r are the top and bottom radiuses
void DrawLtCylinder1(int n, double topr, double bottomr)
{
	double alpha, teta = 2 * PI / n;

	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		glBegin(GL_POLYGON);
		glNormal3d(sin(alpha), bottomr * (bottomr - topr), cos(alpha));
		glVertex3d(topr * sin(alpha), 1, topr * cos(alpha)); //1
		glNormal3d(sin(alpha + teta), bottomr * (bottomr - topr), cos(alpha + teta));

		glVertex3d(topr * sin(alpha + teta), 1, topr * cos(alpha + teta));//2
		glVertex3d(bottomr * sin(alpha + teta), 0, bottomr * cos(alpha + teta));//3
		glNormal3d(sin(alpha), bottomr * (bottomr - topr), cos(alpha));
		glVertex3d(bottomr * sin(alpha), 0, bottomr * cos(alpha));// 4
		glEnd();
	}
}
void DrawCylinder(int n)
{
	double alpha, teta = 2 * PI / n;

	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		glBegin(GL_POLYGON);
		glColor3d(fabs(sin(alpha)) / 2, (1 + cos(alpha)) / 3, fabs(cos(alpha + PI / 3)) / 2);
		glVertex3d(sin(alpha), 1, cos(alpha)); //1
		glColor3d(fabs(sin(alpha + teta)) / 2, (1 + cos(alpha + teta)) / 3, fabs(cos(alpha + teta + PI / 3)) / 2);
		glVertex3d(sin(alpha + teta), 1, cos(alpha + teta));//2
		glColor3d(fabs(sin(alpha)), (1 + cos(alpha)), fabs(cos(alpha)));
		glVertex3d(sin(alpha + teta), 0, cos(alpha + teta));//3
		glVertex3d(sin(alpha), 0, cos(alpha));// 4
		glEnd();
	}
}
void DrawConus(int n)
{
	double alpha, teta = 2 * PI / n;

	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		glBegin(GL_POLYGON);
		glColor3d(fabs(sin(alpha)) / 2, (1 + cos(alpha)) / 3, fabs(cos(alpha + PI / 3)) / 2);
		glVertex3d(0, 1, 0); //1
		glColor3d(fabs(sin(alpha)), (1 + cos(alpha)), fabs(cos(alpha)));
		glVertex3d(sin(alpha + teta), 0, cos(alpha + teta));//3
		glVertex3d(sin(alpha), 0, cos(alpha));// 4
		glEnd();
	}
}
void StrangeCylinder(int n)
{
	double alpha, teta = 2 * PI / n;

	for (alpha = 0; alpha <= 2 * PI; alpha += 2 * teta)
	{
		glBegin(GL_POLYGON);
		glColor3d(fabs(sin(alpha)) / 2, (1 + cos(alpha)) / 3, fabs(cos(alpha + PI / 3)) / 2);
		glVertex3d(sin(alpha), 1, cos(alpha)); //1
		glColor3d(fabs(sin(alpha + teta)) / 2, (1 + cos(alpha + teta)) / 3, fabs(cos(alpha + teta + PI / 3)) / 2);
		glVertex3d(sin(alpha + teta), 1, cos(alpha + teta));//2
		glColor3d(fabs(sin(alpha)), (1 + cos(alpha)), fabs(cos(alpha)));
		glVertex3d(sin(alpha + 9 * teta), 0, cos(alpha + 9 * teta));//3
		glVertex3d(sin(alpha + 8 * teta), 0, cos(alpha + 8 * teta));// 4
		glEnd();
	}
}

void DrawLtSphere(int slices, int stacks)
{
	double beta = PI / slices, gamma; // gamma is running angle from -90 to 90
	int i;
	for (gamma = -PI / 2, i = 0; gamma <= PI / 2; gamma += beta, i++)
	{
		glPushMatrix();
		//		glRotated(i * offset/10, 0, 1, 0);
		glTranslated(0, sin(gamma), 0);
		glScaled(1, sin(gamma + beta) - sin(gamma), 1);
		DrawLtCylinder1(stacks, cos(gamma + beta), cos(gamma));
		glPopMatrix();
	}
}
void HandleShooting() {
	if (shouldShoot && !shouldExplode) {
		glPushMatrix();
		glColor3d(0, 0, 0);
		glTranslated(xf, yf, zf);

		glRotated(180, 0, 1, 0);

		glRotated(yawi * 180 / PI, 0, 1, 0);

		glTranslated(0, 0, -5);
		glTranslated(2.5, 2, -12);
		glutSolidSphere(1, 20, 20);
		glPopMatrix();
	}
	if (shouldExplode) {
		glPushMatrix();
		//glColor3d(0.9, cos(tExplostion*100), sin(tExplostion*100));
		glTranslated(xf, yf, zf);
		glRotated(180, 0, 1, 0);
		glRotated(yawi * 180 / PI, 0, 1, 0);
		glTranslated(0, 0, -5);
		glTranslated(2.5, 2, -12);
		glScaled(2 + tExplostion, 2 + tExplostion, 2 + tExplostion);
		DrawTexSphere(20,20, 7, 1, 1 );
		//glutSolidSphere(2 + tExplostion, 20, 20);
		glPopMatrix();
	}
}
void DrawTank()
{
	glPushMatrix();
	glPushMatrix();
	glRotated(90, 0, 0, 1);
	glRotated(90, 0, 0, -1); // put it to -x axis
	glScaled(1.2, 1, 1.2);

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 9); // new texture name
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // GL_MODULATE

	/////////////////////////////////////////////////////////////////chain left 
	// top face
	glBegin(GL_POLYGON);
	glColor3d(1, 0, 0); // red
	glTexCoord2d(0, 0); glVertex3d(3, 0, 4);
	glTexCoord2d(0, 1); glVertex3d(3, 0, -4);
	glTexCoord2d(1, 1); glVertex3d(1, 0, -4);
	glTexCoord2d(1, 0); glVertex3d(1, 0, 4);
	glEnd();
	// bottom face
	glBegin(GL_POLYGON);
	glColor3d(1, 0, 1); // magenta
	glTexCoord2d(0, 0); glVertex3d(3, -3, 3);
	glTexCoord2d(0, 1); glVertex3d(3, -3, -3);
	glTexCoord2d(1, 1); glVertex3d(1, -3, -3);
	glTexCoord2d(1, 0); glVertex3d(1, -3, 3);
	glEnd();
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 6); // new texture name
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // GL_MODULATE
	// left face 
	glBegin(GL_POLYGON);
	glColor3d(1, 1, 1); // white
	glTexCoord2d(0, 0); glVertex3d(1, 0, 4);
	glTexCoord2d(0, 1); glVertex3d(1, 0, -4);
	glTexCoord2d(1, 1); glVertex3d(1, -3, -3);
	glTexCoord2d(1, 0); glVertex3d(1, -3, 3);
	glEnd();
	// right face
	glBegin(GL_POLYGON);
	glColor3d(1, 1, 0); // yellow
	glTexCoord2d(0, 0); glVertex3d(3, 0, 4);
	glTexCoord2d(0, 1); glVertex3d(3, 0, -4);
	glTexCoord2d(1, 1); glVertex3d(3, -3, -3);
	glTexCoord2d(1, 0); glVertex3d(3, -3, 3);
	glEnd();
	
	
	// back
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 8); // new texture name
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // GL_MODULATE
	glBegin(GL_POLYGON);
	glColor3d(0, 1, 0); // green
	glTexCoord2d(0, 0); glVertex3d(3, 0, -4);
	glTexCoord2d(0, 1); glVertex3d(1, 0, -4);
	glTexCoord2d(1, 1); glVertex3d(1, -3, -3);
	glTexCoord2d(1, 0); glVertex3d(3, -3, -3);
	glEnd();
	// front face
	glBegin(GL_POLYGON);
	glColor3d(0, 0, 4); // blue
	glTexCoord2d(0, 0); glVertex3d(1, 0, 4);
	glTexCoord2d(0, 1); glVertex3d(3, 0, 4);
	glTexCoord2d(1, 1); glVertex3d(3, -3, 3);
	glTexCoord2d(1, 0); glVertex3d(1, -3, 3);
	glEnd();


	/////////////////////////////////////////////////////////////////chain right
		// back
	glBegin(GL_POLYGON);
	glColor3d(0, 1, 0); // green
	glTexCoord2d(0, 0); glVertex3d(-1, 0, -4);
	glTexCoord2d(0, 1); glVertex3d(-3, 0, -4);
	glTexCoord2d(1, 1); glVertex3d(-3, -3, -3);
	glTexCoord2d(1, 0); glVertex3d(-1, -3, -3);
	glEnd();

	// front face
	glBegin(GL_POLYGON);
	glColor3d(0, 0, 4); // blue
	glTexCoord2d(0, 0); glVertex3d(-3, 0, 4);
	glTexCoord2d(0, 1); glVertex3d(-1, 0, 4);
	glTexCoord2d(1, 1); glVertex3d(-1, -3, 3);
	glTexCoord2d(1, 0); glVertex3d(-3, -3, 3);
	glEnd();

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 9); // new texture name
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // GL_MODULATE
// top face
	glBegin(GL_POLYGON);
	glColor3d(1, 0, 0); // red
	glTexCoord2d(0, 0); glVertex3d(-1, 0, 4);
	glTexCoord2d(0, 1); glVertex3d(-1, 0, -4);
	glTexCoord2d(1, 1); glVertex3d(-3, 0, -4);
	glTexCoord2d(1, 0); glVertex3d(-3, 0, 4);
	glEnd();
	// bottom face
	glBegin(GL_POLYGON);
	glColor3d(1, 0, 1); // magenta
	glTexCoord2d(0, 0); glVertex3d(-1, -3, 3);
	glTexCoord2d(0, 1); glVertex3d(-1, -3, -3);
	glTexCoord2d(1, 1); glVertex3d(-3, -3, -3);
	glTexCoord2d(1, 0); glVertex3d(-3, -3, 3);
	glEnd();
	// left face 
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 6); // new texture name
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // GL_MODULATE
	glBegin(GL_POLYGON);
	glColor3d(1, 1, 1); // white
	glTexCoord2d(0, 0); glVertex3d(-3, 0, 4);
	glTexCoord2d(0, 1); glVertex3d(-3, 0, -4);
	glTexCoord2d(1, 1); glVertex3d(-3, -3, -3);
	glTexCoord2d(1, 0); glVertex3d(-3, -3, 3);
	glEnd();
	// right face
	glBegin(GL_POLYGON);
	glColor3d(1, 1, 0); // yellow
	glTexCoord2d(0, 0); glVertex3d(-1, 1, 4);
	glTexCoord2d(0, 1); glVertex3d(-1, 1, -4);
	glTexCoord2d(1, 1); glVertex3d(-1, -3, -3);
	glTexCoord2d(1, 0); glVertex3d(-1, -3, 3);
	glEnd();


	///////////////////////////////////////////////

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 3); // new texture name
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // GL_MODULATE

	// top face
	glBegin(GL_POLYGON);
	//glColor3d(1, 0, 0); // red
	glTexCoord2d(0, 0);	glVertex3d(2, 1, 4);
	glTexCoord2d(0, 1);	glVertex3d(2, 1, -4);
	glTexCoord2d(1, 1);	glVertex3d(-2, 1, -4);
	glTexCoord2d(1, 0);	glVertex3d(-2, 1, 4);
	glEnd();
	// bottom face
	glBegin(GL_POLYGON);
	//glColor3d(1, 0, 1); // magenta
	glTexCoord2d(0, 0);	glVertex3d(2, -1, 4);
	glTexCoord2d(0, 1);	glVertex3d(2, -1, -4);
	glTexCoord2d(1, 1);	glVertex3d(-2, -1, -4);
	glTexCoord2d(1, 0);	glVertex3d(-2, -1, 4);
	glEnd();
	// left face 
	glBegin(GL_POLYGON);
	//glColor3d(1, 1, 1); // white
	glTexCoord2d(0, 0);	glVertex3d(-2, 1, 4);
	glTexCoord2d(0, 1);	glVertex3d(-2, 1, -4);
	glTexCoord2d(1, 1);	glVertex3d(-2, -1, -4);
	glTexCoord2d(1, 0);	glVertex3d(-2, -1, 4);
	glEnd();
	// right face
	glBegin(GL_POLYGON);
	//glColor3d(1, 1, 0); // yellow
	glTexCoord2d(0, 0);	glVertex3d(2, 1, 4);
	glTexCoord2d(0, 1);	glVertex3d(2, 1, -4);
	glTexCoord2d(1, 1);	glVertex3d(2, -1, -4);
	glTexCoord2d(1, 0);	glVertex3d(2, -1, 4);
	glEnd();
	// back
	glBegin(GL_POLYGON);
	//glColor3d(0, 1, 0); // green
	glTexCoord2d(0, 0);	glVertex3d(2, 1, -4);
	glTexCoord2d(0, 1);	glVertex3d(-2, 1, -4);
	glTexCoord2d(1, 1);	glVertex3d(-2, -1, -4);
	glTexCoord2d(1, 0);	glVertex3d(2, -1, -4);
	glEnd();

	// front face
	glBegin(GL_POLYGON);
	//glColor3d(0, 0, 4); // blue
	glTexCoord2d(0, 0);	glVertex3d(-2, 1, 4);
	glTexCoord2d(0, 1);	glVertex3d(2, 1, 4);
	glTexCoord2d(1, 1);	glVertex3d(2, -1, 4);
	glTexCoord2d(1, 0);	glVertex3d(-2, -1, 4);
	glEnd();
	glPopMatrix();
	//Draw top cube
	glPushMatrix();
	glRotated(-apitch2h * 180 * 0.65 / PI, 0, 1, 0);
	glPushMatrix();
	glBegin(GL_POLYGON);
	//glColor3d(1, 0, 1); // magenta
	glTexCoord2d(0, 0);	glVertex3d(1.5, 2, 2);
	glTexCoord2d(0, 1);	glVertex3d(1.5, 2, -2);
	glTexCoord2d(1, 1);	glVertex3d(-1.5, 2, -2);
	glTexCoord2d(1, 0);	glVertex3d(-1.5, 2, 2);
	glEnd();
	// bottom face
	glBegin(GL_POLYGON);
	glColor3d(1, 0, 0); // red
	glTexCoord2d(0, 0);	glVertex3d(1.5, 1, 2);
	glTexCoord2d(0, 1);	glVertex3d(1.5, 1, -2);
	glTexCoord2d(1, 1);	glVertex3d(-1.5, 1, -2);
	glTexCoord2d(1, 0);	glVertex3d(-1.5, 1, 2);
	glEnd();
	// left face 
	glBegin(GL_POLYGON);
	glColor3d(1, 1, 0); // yellow
	glTexCoord2d(0, 0);	glVertex3d(-1.5, 2, 2);
	glTexCoord2d(0, 1);	glVertex3d(-1.5, 2, -2);
	glTexCoord2d(1, 1);	glVertex3d(-1.5, 1, -2);
	glTexCoord2d(1, 0);	glVertex3d(-1.5, 1, 2);
	glEnd();
	// right face
	glBegin(GL_POLYGON);
	glColor3d(1, 1, 1); // white
	glTexCoord2d(0, 0);	glVertex3d(1.5, 2, 2);
	glTexCoord2d(0, 1);	glVertex3d(1.5, 2, -2);
	glTexCoord2d(1, 1);	glVertex3d(1.5, 1, -2);
	glTexCoord2d(1, 0);	glVertex3d(1.5, 1, 2);
	glEnd();
	// back	

	glBegin(GL_POLYGON);
	glColor3d(0, 0, 4); // blue
	glTexCoord2d(0, 0);	glVertex3d(1.5, 2, -2);
	glTexCoord2d(0, 1);	glVertex3d(-1.5, 2, -2);
	glTexCoord2d(1, 1);	glVertex3d(-1.5, 1, -2);
	glTexCoord2d(1, 0);	glVertex3d(1.5, 1, -2);
	glEnd();

	// front face
	glBegin(GL_POLYGON);
	glColor3d(0, 1, 0); // green
	glTexCoord2d(0, 0);	glVertex3d(-1.5, 2, 2);
	glTexCoord2d(0, 1);	glVertex3d(1.5, 2, 2);
	glTexCoord2d(1, 1);	glVertex3d(1.5, 1, 2);
	glTexCoord2d(1, 0);	glVertex3d(-1.5, 1, 2);
	glEnd();

	glPopMatrix();
	//Draw totah

	glPushMatrix();
	glPushMatrix();
	glRotated(-apitch1v * 180 * 0.6 / PI, 1, 0, 0);

	glRotated(90, 1, 0.5, 0.5);
	glTranslated(0, 0, 0);
	glScaled(0.2, 9, 0.2);
	DrawTexCylinder1(400, 1.5, 4, 6, 80);
	glPopMatrix();
	glPopMatrix();

	//End totah
	glPopMatrix();
	// Edd top cube
	glDisable(GL_TEXTURE_2D);
	glPopMatrix();




	
}
void DrawTree2(int x, int z)
{
	glPushMatrix();
	glTranslated(x, 0, z);
	glScaled(1, 4, 1);
	glTranslated(0, 2, 0);
	DrawTexSphere(140, 40, 4, 1, 1);
	glPopMatrix();

	glPushMatrix();
	glTranslated(x, 0, z);
	glScaled(0.1, 8, 0.1);
	DrawTexCylinder1(400, 3, 4, 5, 80);
	glPopMatrix();

}
void DrawPitchControl() // 2D
{
	// background
	glColor3d(0.4, 0.2, 0);
	glBegin(GL_POLYGON);
	glVertex2d(-1, -1);
	glVertex2d(-1, 1);
	glVertex2d(1, 1);
	glVertex2d(1, -1);
	glEnd();

	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(0, -1);
	glVertex2d(0, 1);
	glEnd();

	// slider
	glColor3d(0.5, 0.5, 0.5);
	glBegin(GL_POLYGON);
	glVertex2d(-0.2, apitch1v - 0.2);
	glVertex2d(-0.2, apitch1v + 0.2);
	glVertex2d(0.2, apitch1v + 0.2);
	glVertex2d(0.2, apitch1v - 0.2);
	glEnd();
	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(-0.2, apitch1v);
	glVertex2d(0.2, apitch1v);
	glEnd();

}
void DrawPitchControl2() // 2D
{
	// background
	glColor3d(0.8, 0.2, 0.3);
	glBegin(GL_POLYGON);
	glVertex2d(-1, -1);
	glVertex2d(-1, 1);
	glVertex2d(1, 1);
	glVertex2d(1, -1);
	glEnd();

	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(-1, 0);
	glVertex2d(1, 0);
	glEnd();

	// slider
	glColor3d(0.5, 0.5, 0.5);
	glBegin(GL_POLYGON);
	glVertex2d(apitch2h - 0.2, -0.2);
	glVertex2d(apitch2h + 0.2, -0.2);
	glVertex2d(apitch2h + 0.2, 0.2);
	glVertex2d(apitch2h - 0.2, 0.2);
	glEnd();
	glColor3d(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2d(apitch2h, -0.2);
	glVertex2d(apitch2h, 0.2);
	glEnd();
}
/*
void DrawCompassControl() // 2D
{
	double alpha, teta = PI/20;
	// background
	glColor3d(0.4, 0.2, 0);
	glBegin(GL_POLYGON);
	glVertex2d(-1, -1);
	glVertex2d(-1, 1);
	glVertex2d(1, 1);
	glVertex2d(1, -1);
	glEnd();

	glColor3d(1, 1, 0.7);
	glBegin(GL_POLYGON);
	for (alpha = 0; alpha <= 2 * PI; alpha += teta)
	{
		glVertex2d(cos(alpha), sin(alpha));
	}
	glEnd();

	// arrow v.1 with sin/cos
//	glColor3d(0.4, 0.2, 0);
//	glBegin(GL_LINES);
//	glVertex2d(0, 0);
//	glVertex2d(sin(yaw),cos(yaw) );
//	glEnd();

// arrow v.2  with rotation
	glColor3d(0.4, 0.2, 0);
	glPushMatrix();
	glRotated(yaw * 180 / PI, 0, 0, -1);
	glBegin(GL_LINE_STRIP);
	glVertex2d(0, 0);
	glVertex2d(0, 1);
	glVertex2d(-0.2, 0.8);
	glVertex2d(0.2, 0.8);
	glVertex2d(0, 1);
	glEnd();
	glPopMatrix();
}
*/
void DrawTree1(int x, int z)
{
	glPushMatrix();
	glTranslated(x, 0, z);
	glScaled(2, 2, 2);
	glTranslated(0, 4, 0);
	DrawTexSphere(140, 40, 4, 1, 1);
	glPopMatrix();
	glPushMatrix();
	glTranslated(x, 0, z);
	glScaled(0.2, 8, 0.2);
	DrawTexCylinder1(400, 1.5, 4, 5, 80);
	glPopMatrix();
}
void SetTree() {

	for (int i = 0; i < 20; i++) {
		DrawTree1(Trees[i].x,Trees[i].z);
	}
	for (int j = 20; j < 40; j++) {
		DrawTree2(Trees[j].x, Trees[j].z);
	}
	
}

// redraw function
void display()
{
	// clean frame buffer and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, W, H);
	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(eyex, eyey, eyez,  // eye
		eyex + dir[0], eyey + sin(pitch), eyez + dir[2],// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	DrawAxes();
	DrawFloor();
	SetTree();
	HandleShooting(); 

	glPushMatrix();
	//while (ground[][]<=2) {
	glTranslated(ax, ay, az);
	//}
	//glTranslated(ax, ay - 9, az);
	glRotated(yaw * 180 / PI, 0, 1, 0);//yaw is in radians so we transform it to degrees
	glRotated(aangular_speed * 5000, 0, 0, -1); // roll around main axis of airplane
	//glRotated(-apitch * 180 / PI, 1, 0, 0);
	DrawTank();
	glPopMatrix();
	glPushMatrix();
	glRotated(offset / 20, 0, 1, 0);
	glTranslated(0, 30, 0);
	glScaled(100, 100, 100);

	//sky
	DrawTexSphere(200, 200, 2, 1, 1);
	glPopMatrix();

	// add pitch control (2D)
	glViewport(W - 100, 0, 100, 100);
	glDisable(GL_DEPTH_TEST); // last painted occludes preveous painted (2D)

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity(); // start transformations from the start
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//Left control
	DrawPitchControl();

	// compass control
	glViewport(0, 0, 100, 100);

	// Right control
	DrawPitchControl2();

	glEnable(GL_DEPTH_TEST); // near occludes far (3D)

	/// <summary>
	/// explotion

	/// </summary>

	glutSwapBuffers(); // show all


}
void combined_display()
{
	// clean frame buffer and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// regular view
	glViewport(0, H / 2, W / 2, H / 2);

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(eyex, eyey, eyez,  // eye
		eyex + dir[0], eyey + sin(pitch), eyez + dir[2],// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	DrawFloor();
	DrawAxes();
	SetTree();
	HandleShooting();

	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(yaw * 180 / PI, 0, 1, 0);//yaw is in radians so we transform it to degrees
	glRotated(aangular_speed * 5000, 0, 0, -1); // roll around main axis of airplane
	glRotated(-apitch1v * 180 / PI, 1, 0, 0);
	DrawTank();

	glPopMatrix();

	// top view
	glViewport(W / 2, H / 2, W / 2, H / 2);
	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(eyex, eyey + 40, eyez,  // eye
		eyex, eyey, eyez - 0.1,// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	DrawFloor();
	DrawAxes();
	SetTree();
	HandleShooting();

	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(yaw * 180 / PI, 0, 1, 0);//yaw is in radians so we transform it to degrees
	glRotated(aangular_speed * 5000, 0, 0, -1); // roll around main axis of airplane
	glRotated(-apitch1v * 180 / PI, 1, 0, 0);
	DrawTank();
	glPopMatrix();

	// cocpit view
	glViewport(0, 0, W, H / 2);
	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(ax + 3 * adir[0], ay + 2 + 3 * adir[1], az + 3 * adir[2],  // eye is now placed into airplane
		ax + 4 * adir[0], ay + 2 + 4 * adir[1], az + 4 * adir[2],// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	DrawFloor();
	DrawAxes();
	SetTree();
	HandleShooting();

	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(yaw * 180 / PI, 0, 1, 0);//yaw is in radians so we transform it to degrees
	glRotated(aangular_speed * 5000, 0, 0, -1); // roll around main axis of airplane
	glRotated(-apitch1v * 180 / PI, 1, 0, 0);
	DrawTank();
	glPopMatrix();

	// add pitch control (2D)
	glViewport(W - 100, 0, 100, 100);
	glDisable(GL_DEPTH_TEST); // last painted occludes preveous painted (2D)

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity(); // start transformations from the start
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//Left control
	DrawPitchControl();

	// compass control
	glViewport(0, 0, 100, 100);

	// Right control
	DrawPitchControl2();

	glEnable(GL_DEPTH_TEST); // near occludes far (3D)
	glutSwapBuffers(); // show all
}
void top_display()
{
	// clean frame buffer and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, W, H);

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(eyex, eyey + 40, eyez,  // eye
		eyex, eyey, eyez - 0.1,// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	DrawFloor();
	DrawAxes();
	SetTree();
	HandleShooting();

	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(yaw * 180 / PI, 0, 1, 0);//yaw is in radians so we transform it to degrees
	glRotated(aangular_speed * 5000, 0, 0, -1); // roll around main axis of airplane
	DrawTank();


	glPopMatrix();
	glutSwapBuffers(); // show all
}
void cockpit_display()
{
	// clean frame buffer and Z-buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, W, H);

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 0.7, 300);
	gluLookAt(ax + 3 * adir[0], ay + 2 + 3 * adir[1], az + 3 * adir[2],  // eye is now placed into airplane
		ax + 4 * adir[0], ay + 2 + 4 * adir[1], az + 4 * adir[2],// center or PointOfInterest
		0, 1, 0); // up

	glMatrixMode(GL_MODELVIEW); // setting the transformation matrix
	glLoadIdentity(); // start transformations fro identity matrix

	DrawFloor();
	DrawAxes();
	SetTree();
	HandleShooting();

	glPushMatrix();
	glTranslated(ax, ay, az);
	glRotated(yaw * 180 / PI, 0, 1, 0);//yaw is in radians so we transform it to degrees
	glRotated(aangular_speed * 5000, 0, 0, -1); // roll around main axis of airplane
	//glRotated(-apitch1v * 180 / PI, 1, 0, 0);
	DrawTank();
	glPopMatrix();

	// add pitch control (2D)
	glViewport(W - 100, 0, 100, 100);
	glDisable(GL_DEPTH_TEST); // last painted occludes preveous painted (2D)

	glMatrixMode(GL_PROJECTION); // setting the projection matrix
	glLoadIdentity(); // start transformations from the start
	glOrtho(-1, 1, -1, 1, -1, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//Left control
	DrawPitchControl();

	// compass control
	glViewport(0, 0, 100, 100);

	// Right control
	DrawPitchControl2();

	glEnable(GL_DEPTH_TEST); // near occludes far (3D)
	glutSwapBuffers(); // show all
}

// runs when nothing happens
void idle()
{
	t += 0.01;
	if (shouldShoot) {
		// calculates distance for x,y,z positions after t seconds of shot
		xf = xi + sin(yawi) * t * vx - horazintalApitch * a * t * t;
		yf = yi + (vy * t * sin(verticalApitch)) + (a * t * t);
		zf = zi + cos(yawi) * t * vz;

		int i, j;

		// iterate over height terrain pixels and check if shot hit terrain height
		for (i = 0; i < GSZ; i++) {
			bool collider = false;
			for (j = 0; j < GSZ; j++) {
				if (ground[i][j] >= yf) {
					shouldShoot = false;
					collider = true;
					break;
				}
			}

			// collision occured, trigger explosion
			if (collider) {
				shouldExplode = true;
				tExplostion = 0;
			}
		}
	}
	// explosion timer
	if (shouldExplode) {
		tExplostion += 0.2;
		if (tExplostion >= 10) {
			shouldExplode = false;
		}
	}

	offset += 0.5;
	
	// ego motion
	sight_angle += angular_speed;

	dir[0] = sin(sight_angle); // x
	dir[1] = sin(pitch);   //y
	dir[2] = cos(sight_angle); // z

	eyex += speed * dir[0];
	eyey += speed * dir[1];
	eyez += speed * dir[2];

	// tank
	yaw += aangular_speed;
	adir[0] = sin(yaw);
	adir[1] = sin(apitch1v);
	adir[2] = cos(yaw);

	ax += aspeed * adir[0];
	ay += aspeed * adir[1];
	az += aspeed * adir[2];

	glutPostRedisplay(); // indirect call to display
}

void SpecialKey(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_UP:

		speed += 0.001;
		break;
	case GLUT_KEY_DOWN:
		speed -= 0.001;
		break;
	case GLUT_KEY_LEFT:
		angular_speed += 0.001;
		break;
	case GLUT_KEY_RIGHT:
		angular_speed -= 0.001;
		break;
	case GLUT_KEY_PAGE_UP:
		pitch += 0.1;
		break;
	case GLUT_KEY_PAGE_DOWN:
		pitch -= 0.1;
		break;
	}
}

void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		//	is_building = !is_building;
		Smooth();
	}
}

void keyboard(unsigned char key, int x, int y)
{

	switch (key)
	{
	case 'a':
		aangular_speed += 0.001;
		break;
	case 'd':
		aangular_speed -= 0.001;
		break;
	case 'w':
		aspeed += 0.001;
		break;
	case 's':
		aspeed -= 0.001;
		break;
	case ' ':
		if (!shouldExplode && !shouldShoot) {
			printf("sdnfjksdnfsd");
			t = 0;
			shouldShoot = true;
			yawi = yaw;
			verticalApitch = apitch1v;
			horazintalApitch = apitch2h;
			xi = ax + 10 * horazintalApitch;
			yi = ay + 10 * verticalApitch;
			zi = az;
		}
		break;
	}
}

void menu(int choice)
{
	int cx = -30 + rand() % 30;
	int cz = -30 + rand() % 30;
	switch (choice)
	{
	case 1:// regular view
		glutDisplayFunc(display);
		break;
	case 2:// top view
		glutDisplayFunc(top_display);
		break;
	case 3:// cockpit view
		glutDisplayFunc(cockpit_display);
		break;
	case 4:// combined view
		glutDisplayFunc(combined_display);
		break;
	}
}

void mouse_motion(int x, int y)
{
	//Right control
	double x1, y1;

	x1 = (2.0 * x / W) - 1; // xx is in range (-1,1)
	y1 = (2.0 * (H - y) / 100) - 1; // yy is in range (-1,1)

	if (0.8 < x1 && x1 < 0.86 && 0 < y1 && y1 < 0.8)
		if (fabs(y1) < 1)
		{
			apitch1v = y1;
		}

	//Left control
	double x2, y2;

	x2 = (2.0 * x / 100) - 1; // xx is in range (-1,1)
	y2 = (2.0 * (H - y) / 100) - 1; // yy is in range (-1,1)

	if (-0.8 < x2 && x2 < 0.8 && apitch2h - 0.8 < y2 && y2 < apitch2h + 0.8)
		if (fabs(x2) < 1) {
			apitch2h = x2;
		}


}

void main(int argc, char* argv[])
{
	glutInit(&argc, argv); // 
	// defines matrices: 1. Color matrix (frame buffer), 
	// 2. video buffer
	// 3. Z-buffer
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(W, H); // physical size of window (in pixels)
	glutInitWindowPosition(200, 100);
	glutCreateWindow("First Example");

	glutDisplayFunc(display); //  display is a refresh window function
	glutIdleFunc(idle); // kind of timer function
//	glutMouseFunc(mouse);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(SpecialKey);
	glutMotionFunc(mouse_motion);
	// menu
	glutCreateMenu(menu);
	glutAddMenuEntry("Regular View", 1);
	glutAddMenuEntry("Top View", 2);
	glutAddMenuEntry("Cockpit View", 3);
	glutAddMenuEntry("Combined View", 4);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	init();
	glutMainLoop();
}