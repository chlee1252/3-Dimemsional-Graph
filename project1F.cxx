/*
Author: Marc Lee
Description: Make 1000 different loaction of 3D objects
*/

#include <iostream>
#include <sstream>                  //used to set different file name
#include <iomanip>                  //used to set different file name
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

#define NORMALS

using std::cerr;
using std::endl;

//int triangleID = -1; //debug code

double ceil_441(double f)
{
    return ceil(f-0.00001);
}

double floor_441(double f)
{
    return floor(f+0.00001);
}


vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };
  

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

double lerp(double v, double v1, double v2, double A, double B)
{    
    return A + ((v-v1)/(v2-v1)) * (B-A);
}

double dot_product(double* A, double* B)
{
    //Lecture 6, pg 32
    return (A[0] * B[0]) + (A[1] * B[1]) + (A[2] * B[2]);
}

double* cross_product(double* A, double* B)
{
    //Lecture 6, pg 33 //0 = x, 1 = y, 2 = z
    double* cp = new double[3];

    cp[0] = A[1] * B[2] - A[2] * B[1];
    cp[1] = A[2] * B[0] - A[0] * B[2];
    cp[2] = A[0] * B[1] - A[1] * B[0];

    return cp;
}

double* normalize(double x, double y, double z)
{
    // ||A|| = sqrt(x^2 + y^2 + z^2)
    //Normal = A / ||A||
    double* normal = new double[3];
    double length = sqrt( (x*x + y*y + z*z));

    normal[0] = x / length; //x
    normal[1] = y / length; //y
    normal[2] = z / length; //z

    return normal;
}

double CalculatePhongShading(LightingParameters &, double *viewDirection, double *normal)
{
    double* light = normalize(lp.lightDir[0], lp.lightDir[1], lp.lightDir[2]);

    double ambient = 1;

    double dot_diffuse = dot_product(light, normal);
    double diffuse = fabs(dot_diffuse);
    //double diffuse = 0;

    double* R = new double[3];
    for (int i = 0; i < 3; i++)
    {
        R[i] = 2 * dot_diffuse * normal[i] - light[i];
    }
    R = normalize(R[0], R[1], R[2]);
    double* V = normalize(viewDirection[0], viewDirection[1], viewDirection[2]);

    double specular = fmax(0, pow(dot_product(R, V), lp.alpha));
    //double specular = 0;
    return (lp.Ka * ambient) + (lp.Kd * diffuse) + (lp.Ks * specular);
}

class Matrix
{
  public:
    double          A[4][4];  // A[i][j] means row i, column j

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void
Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix
Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void
Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Screen
{
  public:
      unsigned char   *buffer;
      double          *zbuffer;
      int width, height;

      void setzBuffer();
      void ImageColor(double r, double c, double colors[3], double shading);

};

void Screen::setzBuffer()
{
    zbuffer = new double[width * height];

    for (int i = 0; i < width * height; i++)
    {
        zbuffer[i] = -1;
    }
}

void Screen::ImageColor(double r, double c, double colors[3], double shading)
{
    for (int i = 0; i < 3; i++)
    {
        colors[i] = colors[i] * shading;
        colors[i] = fmin(colors[i], 1);
    }

    int index = r * width + c;
    buffer[3*index+0] = ceil_441(colors[0] * 255);
    buffer[3*index+1] = ceil_441(colors[1] * 255);
    buffer[3*index+2] = ceil_441(colors[2] * 255);

}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];

    double*          u;
    double*          v;
    double*          w;

    void            positionSetup();

    Matrix          ViewTransform(void);
    Matrix          CameraTransform(void);
    Matrix          DeviceTransform(Screen s);
};

void Camera::positionSetup()
{
    //Lecture 6, pg 46
    double* new_up = new double[3];
    u = new double[3];
    v = new double[3];
    w = new double[3];

    for(int i = 0; i < 3; i++)
    {
        w[i] = position[i] - focus[i];
        new_up[i] = up[i];
    }

    u = cross_product(new_up, w);
    v = cross_product(w, u);

    u = normalize(u[0], u[1], u[2]);
    v = normalize(v[0], v[1], v[2]);
    w = normalize(w[0], w[1], w[2]);

    delete new_up;

}

Matrix Camera::ViewTransform(void)
{
    //Lecture 6, pg 61
    Matrix Transform_view;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Transform_view.A[i][j] = 0;
        }
    }

    Transform_view.A[0][0] = 1 / (tan(angle / 2));

    Transform_view.A[1][1] = 1 / (tan(angle / 2));

    Transform_view.A[2][2] = (far + near) / (far - near);
    Transform_view.A[2][3] = -1;

    Transform_view.A[3][2] = (2 * far * near) / (far - near);

    return Transform_view;
}

Matrix Camera::CameraTransform(void)
{
    //Lecture 6, pg 56
    Matrix Transform_camera;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Transform_camera.A[i][j] = 0;
        }
    }

    //t = (0,0,0) - position
    double* t = new double[3];
    t[0] = 0 - position[0];
    t[1] = 0 - position[1];
    t[2] = 0 - position[2];

    Transform_camera.A[0][0] = u[0];  //u.x
    Transform_camera.A[0][1] = v[0];  //v.x
    Transform_camera.A[0][2] = w[0];  //w.x

    Transform_camera.A[1][0] = u[1];  //u.y
    Transform_camera.A[1][1] = v[1];  //v.y
    Transform_camera.A[1][2] = w[1];  //w.y

    Transform_camera.A[2][0] = u[2];  //u.z
    Transform_camera.A[2][1] = v[2];  //v.z
    Transform_camera.A[2][2] = w[2];  //w.z

    Transform_camera.A[3][0] = dot_product(u, t);  //u dot t
    Transform_camera.A[3][1] = dot_product(v, t);  //v dot t
    Transform_camera.A[3][2] = dot_product(w, t);  //w dot t
    Transform_camera.A[3][3] = 1;

    delete t;
    return Transform_camera;
}

Matrix Camera::DeviceTransform(Screen s)
{
    //Lecture 6, pg 28
    Matrix Transform_Device;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Transform_Device.A[i][j] = 0;
        }
    }

    double x = s.width / 2;
    double y = s.height / 2;
    Transform_Device.A[0][0] = x;

    Transform_Device.A[1][1] = y;

    Transform_Device.A[2][2] = 1;

    Transform_Device.A[3][0] = x;
    Transform_Device.A[3][1] = y;
    Transform_Device.A[3][3] = 1;

    return Transform_Device;
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      double         colors[3][3];
      double         normals[3][3];
      double         shading[3];

      int right, left, min, mid, max;
      double slopeL, slopeR, slopeM;
      double yinterL, yinterR, yinterM;

      double* ptIn;
      double* ptOut;

      void setPoints();

      void setRight_Bottom();
      void setLeft_Bottom();

      void setRight_Top();
      void setLeft_Top();

      double Bottom_leftL(double y);
      double Bottom_rightL(double y);

      double Top_leftL(double y);
      double Top_rightL(double y);

      double MidLine(double y);

      void TransformTrianglesToDeviceSpace(Matrix x);

      void CalculateShading(Camera c);

};

void Triangle::setLeft_Bottom()
{
    if (X[mid] < X[max])
    {
        left = mid;
    }
    else
    {
        left = max;
    }
}

void Triangle::setRight_Bottom()
{
    if (X[mid] < X[max])
    {
        right = max;
    }
    else
    {
        right = mid;
    }
}

void Triangle::setLeft_Top()
{
    if (X[min] < X[mid])
    {
        left = min;
    }
    else
    {
        left = mid;
    }

}

void Triangle::setRight_Top()
{
    if (X[min] < X[mid])
    {
        right = mid;
    }
    else
    {
        right = min;
    }

}

void Triangle::setPoints()
{
    if (Y[0] > Y[1])
    {
        if (Y[0] > Y[2])
        {
            max = 0;
            if (Y[1] > Y[2])
            {
                mid = 1;
                min = 2;
            }
            else {
                mid = 2;
                min = 1;
            }
        }
        else
        {
            mid = 0;
            if ( Y[1] > Y[2])
            {
                min = 2;
                max = 1;
            }
            else
            {
                max = 2;
                min = 1;
            }
        }
    }
    else
    {
        if (Y[1] > Y[2])
        {
            max = 1;
            if (Y[0] > Y[2])
            {
                mid = 0;
                min = 2;
            }
            else
            {
                mid = 2;
                min = 0;
            }
        }
        else
        {
            min = 0;
            mid = 1;
            max = 2;
        }
    }
}

double Triangle::Bottom_rightL(double y)
{
    setRight_Bottom();

    double x1 = X[min];
    double y1 = Y[min];
    double x2 = X[right];
    double y2 = Y[right];

    if (x1 == x2)
    {
        slopeR = x1;
        return slopeR;
    }
    else 
    {
        slopeR = (y2 - y1) / (x2 - x1);
        yinterR = y1 - slopeR * x1;
        return (y - yinterR) / slopeR;
    }
}

double Triangle::Bottom_leftL(double y)
{
    setLeft_Bottom();
    double x1 = X[min];
    double y1 = Y[min];
    double x2 = X[left];
    double y2 = Y[left];

    if (x1 == x2)
    {
        slopeL = x1;
        return slopeL;
    }
    else
    {
        slopeL = (y2 - y1) / (x2 - x1);
        yinterL = y1 - slopeL * x1;
        return (y - yinterL) / slopeL;
    }

}

double Triangle::Top_leftL(double y)
{
    setLeft_Top();
    double x1 = X[left];
    double y1 = Y[left];
    double x2 = X[max];
    double y2 = Y[max];

    if (x1 == x2)
    {
        slopeL = x1;
        return slopeL;
    }
    else
    {
        slopeL = (y2 - y1) / (x2 - x1);
        yinterL = y1 - slopeL * x1;
        return (y - yinterL) / slopeL;
    }
}

double Triangle::Top_rightL(double y)
{
    setRight_Top();
    double x1 = X[right];
    double y1 = Y[right];
    double x2 = X[max];
    double y2 = Y[max];

    if (x1 == x2)
    {
        slopeR = x1;
        return slopeR;
    }
    else
    {
        slopeR = (y2 - y1) / (x2 - x1);
        yinterR = y1 - slopeR * x1;
        return (y - yinterR) / slopeR;
    }
}

double Triangle::MidLine(double y)
{
    double x1 = X[min];
    double y1 = Y[min];
    double x2 = X[max];
    double y2 = Y[max];

    if (x1 == x2)
    {
        slopeM = x1;
        return slopeM;
    }
    else 
    {
        slopeM = (y2 - y1) / (x2 - x1);
        yinterM = y1 - slopeM * x1;
        return (y - yinterM) / slopeM;
    }
}

void Triangle::TransformTrianglesToDeviceSpace(Matrix x)
{
    //Lecture 6, pg 37
    for (int i = 0; i < 3; i++)
    {
        ptIn = new double[4];
        ptOut = new double[4];

        ptIn[0] = X[i];
        ptIn[1] = Y[i];
        ptIn[2] = Z[i];
        ptIn[3] = 1;

        x.TransformPoint(ptIn, ptOut);
        
        //divide by w, Lecture 6, pg 39
        ptOut[0] = ptOut[0] / ptOut[3];   
        ptOut[1] = ptOut[1] / ptOut[3];
        ptOut[2] = ptOut[2] / ptOut[3];

        X[i] = ptOut[0];
        Y[i] = ptOut[1];
        Z[i] = ptOut[2];

        delete ptIn;
        delete ptOut;
    }
}

void Triangle::CalculateShading(Camera c)
{
    for (int i = 0; i < 3; i++)
    {
        double* pos = new double[3];
        pos[0] = c.position[0] - X[i];
        pos[1] = c.position[1] - Y[i];
        pos[2] = c.position[2] - Z[i];

        shading[i] = CalculatePhongShading(lp, pos, normals[i]);

        delete pos;
    }
}

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera
GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

double* setColor(double v, double v1, double v2, double* A, double* B)
{
    double* color = new double[3];

    color[0] = lerp(v, v1, v2, A[0], B[0]);
    color[1] = lerp(v, v1, v2, A[1], B[1]);
    color[2] = lerp(v, v1, v2, A[2], B[2]);

    return color;
}

void RasterizeGoingDownTriangle(Triangle triangles, Screen screen)
{
    double rowMin = ceil_441(triangles.Y[triangles.min]);
    double rowMax = floor_441(triangles.Y[triangles.max]);

    for (int r = rowMin; r <= rowMax; r++)
    {
        double leftEnd = triangles.Bottom_leftL(r);
        double rightEnd = triangles.Bottom_rightL(r);

        int right = triangles.right;
        int left = triangles.left;
        int min = triangles.min;

        double zleft = lerp(r, triangles.Y[left], triangles.Y[min], triangles.Z[left], triangles.Z[min]);
        double zright = lerp(r, triangles.Y[right], triangles.Y[min], triangles.Z[right], triangles.Z[min]);

        double sleft = lerp(r, triangles.Y[left], triangles.Y[min], triangles.shading[left], triangles.shading[min]);
        double sright = lerp(r, triangles.Y[right], triangles.Y[min], triangles.shading[right], triangles.shading[min]);

        double* leftC = setColor(r, triangles.Y[left], triangles.Y[min], triangles.colors[left], triangles.colors[min]);
        double* rightC = setColor(r, triangles.Y[right], triangles.Y[min], triangles.colors[right], triangles.colors[min]);

        for (int c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++)
        {
            double* color = setColor(c, leftEnd, rightEnd, leftC, rightC);
            double z = lerp(c, leftEnd, rightEnd, zleft ,zright);
            double shade = lerp(c, leftEnd, rightEnd, sleft, sright);
            if (r >= 0 && r < screen.height && c >= 0 && c < screen.width)
            {              
                if (z > screen.zbuffer[r * screen.width + c])
                {
                    screen.zbuffer[r * screen.width + c] = z;
                    screen.ImageColor(r, c, color, shade);
                }
            }
            delete color;
        }
        delete leftC;
        delete rightC;
    }

}

void RasterizeGoingUpTriangle(Triangle triangles, Screen screen)
{
    double rowMin = ceil_441(triangles.Y[triangles.min]);
    double rowMax = floor_441(triangles.Y[triangles.max]);

    for (int r = rowMin; r <= rowMax; r++)
    {
        double leftEnd = triangles.Top_leftL(r);
        double rightEnd = triangles.Top_rightL(r);

        int right = triangles.right;
        int left = triangles.left;
        int max = triangles.max;

        double zleft = lerp(r, triangles.Y[left], triangles.Y[max], triangles.Z[left], triangles.Z[max]);
        double zright = lerp(r, triangles.Y[right], triangles.Y[max], triangles.Z[right], triangles.Z[max]);

        double sleft = lerp(r, triangles.Y[left], triangles.Y[max], triangles.shading[left], triangles.shading[max]);
        double sright = lerp(r, triangles.Y[right], triangles.Y[max], triangles.shading[right], triangles.shading[max]);

        double* leftC = setColor(r, triangles.Y[left], triangles.Y[max], triangles.colors[left], triangles.colors[max]);
        double* rightC = setColor(r, triangles.Y[right], triangles.Y[max], triangles.colors[right], triangles.colors[max]);

        for (int c = ceil_441(leftEnd); c <= floor_441(rightEnd); c++)
        {
            double* color = setColor(c, leftEnd, rightEnd, leftC, rightC);
            double z = lerp(c, leftEnd, rightEnd, zleft, zright);
            double shade = lerp(c, leftEnd, rightEnd, sleft, sright);
            if (r >= 0 && r < screen.height && c >= 0 && c < screen.width)
            {              
                if (z > screen.zbuffer[r * screen.width + c])
                {
                    screen.zbuffer[r * screen.width + c] = z;
                    screen.ImageColor(r, c, color, shade);
                }
            }
            delete color;
        }
        delete leftC;
        delete rightC;
    }
}

void RasterizeArbitraryTriangle(Triangle triangles, Screen screen)
{
    Triangle topT = triangles;
    Triangle botT = triangles;

    double new_x = triangles.MidLine(triangles.Y[triangles.mid]);
    double new_y = triangles.Y[triangles.mid];
    int min = triangles.min;
    int max = triangles.max;
    
    botT.colors[max][0] = lerp(new_y, botT.Y[min], botT.Y[max], botT.colors[min][0], botT.colors[max][0]);
    botT.colors[max][1] = lerp(new_y, botT.Y[min], botT.Y[max], botT.colors[min][1], botT.colors[max][1]);
    botT.colors[max][2] = lerp(new_y, botT.Y[min], botT.Y[max], botT.colors[min][2], botT.colors[max][2]);
    botT.Z[max] = lerp(new_y, botT.Y[min], botT.Y[max], botT.Z[min], botT.Z[max]);

    botT.shading[max] = lerp(new_y, botT.Y[min], botT.Y[max], botT.shading[min], botT.shading[max]);

    botT.X[botT.max] = new_x;
    botT.Y[botT.max] = new_y;
    //--------------------------------------------------------------------------
    topT.colors[min][0] = lerp(new_y, topT.Y[min], topT.Y[max], topT.colors[min][0], topT.colors[max][0]);
    topT.colors[min][1] = lerp(new_y, topT.Y[min], topT.Y[max], topT.colors[min][1], topT.colors[max][1]);
    topT.colors[min][2] = lerp(new_y, topT.Y[min], topT.Y[max], topT.colors[min][2], topT.colors[max][2]);
    topT.Z[min] = lerp(new_y, topT.Y[min], topT.Y[max], topT.Z[min], topT.Z[max]);

    topT.shading[min] = lerp(new_y, topT.Y[min], topT.Y[max], topT.shading[min], topT.shading[max]);

    topT.X[topT.min] = new_x;
    topT.Y[topT.min] = new_y;

    RasterizeGoingUpTriangle(topT, screen);
    RasterizeGoingDownTriangle(botT, screen);
    

}

bool isTriangle(Triangle triangles)
{
    if (triangles.X[0] == triangles.X[1] && triangles.Y[0] == triangles.Y[1])
        return false;
    if (triangles.X[0] == triangles.X[2] && triangles.Y[0] == triangles.Y[2])
        return false;
    if (triangles.X[1] == triangles.X[2] && triangles.Y[1] == triangles.Y[2])
        return false;
    
    return true;
}

void RasterizeTriangle(Triangle triangles, Screen screen)
{
    triangles.setPoints();
    if (isTriangle(triangles))
    {
        if (triangles.Y[triangles.min] == triangles.Y[triangles.mid])
        {
            RasterizeGoingUpTriangle(triangles, screen);
            return;
        }

        if (triangles.Y[triangles.max] == triangles.Y[triangles.mid])
        {
            RasterizeGoingDownTriangle(triangles, screen);
            return;
        }

        RasterizeArbitraryTriangle(triangles, screen);
        return;
    }
}

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

Screen InitializeScreen(unsigned char *buffer)
{
    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;
    screen.setzBuffer();

    int npixels = 1000 * 1000;
    for (int i = 0; i < npixels*3; i++)
    {
        screen.buffer[i] = 0;
    }

    return screen;
}

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   std::vector<Triangle> triangles = GetTriangles();  
   
   Screen screen;

   for (int i = 0; i < 4; i++)
   {
        int v = 250 * i;
        screen = InitializeScreen(buffer);
        Camera c = GetCamera(i, 1000);
        c.positionSetup();
        Matrix ct = c.CameraTransform();
        Matrix vt = c.ViewTransform();
        Matrix dt = c.DeviceTransform(screen);

        //Camera -> View -> Device
        Matrix compose = Matrix::ComposeMatrices(Matrix::ComposeMatrices(ct, vt), dt);

        for (int j = 0; j < triangles.size(); j++) 
        {
                //triangleID = i;             //debug code
                triangles[j].setPoints();
                triangles[j].CalculateShading(c);
                Triangle t = triangles[j];    //save original triangle data
                triangles[j].TransformTrianglesToDeviceSpace(compose);
                RasterizeTriangle(triangles[j], screen);
                triangles[j] = t;             //restore original triangle data for later use
        }

        if (i == 0)
        {
            WriteImage(image, "frame000");
        }

        std::stringstream ss;
        ss << std::setw(3) << std::setfill('0') << i;
        std::string file = "subImageDir/frame" + ss.str();
        const char *name = file.c_str();
        WriteImage(image, name);
   }

   delete[] screen.zbuffer;
}
