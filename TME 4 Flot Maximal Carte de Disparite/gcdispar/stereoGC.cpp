// Imagine++ project
// Project:  GraphCuts
// Author:   Renaud Marlet

#include <Imagine/Images.h>
#include <iostream>
#include <algorithm>
#include <map>
#include "maxflow/graph.h"

using namespace std;
using namespace Imagine;

typedef Image<byte> byteImage;
typedef Image<double> doubleImage;

// Return image of mean intensity value over (2n+1)x(2n+1) patch
doubleImage meanImage(const doubleImage& I, int n) {
    // Create image for mean values
    int w = I.width(), h = I.height();
    doubleImage IM(w,h);
    // Compute patch area
    double area = (2*n+1)*(2*n+1);
    // For each pixel
    for (int i=0; i<w; i++)
        for (int j=0; j<h; j++) {
            // If pixel is close to border (<n pixels) mean is meaningless
            if (j-n<0 || j+n>=h || i-n<0 ||i+n>=w) {
                IM(i,j)=0;
                continue;
            }
            double sum = 0;
            for(int x = i-n; x <= i+n; x++)
                for (int y = j-n; y <= j+n; y++)
                    sum += I(x,y);
            IM(i,j) = sum / area;
        }
    return IM;
}

// Compute correlation between two pixels in images 1 and 2
double correl(const doubleImage& I1,  // Image 1
              const doubleImage& I1M, // Image of mean value over patch
              const doubleImage& I2,  // Image2
              const doubleImage& I2M, // Image of mean value over patch
              int u1, int v1,         // Pixel of interest in image 1
              int u2, int v2,         // Pixel of interest in image 2
              int n) {                // Half patch size
    // Initialize correlation
    double c = 0;
    // For each pixel displacement in patch
    for (int x=-n; x<=n; x++) {
        for(int y=-n; y<=n; y++){
            if(u1+x<I1.width() && u1+x>=0 && 
                v1+y<I1.height()  && v1+y>=0 &&  
                u2+x<I2.width() && u2+x>=0 && 
                v2+y<I2.height()  && v2+y>=0 && 
                u1<I1.width() && u1>=0 && 
                u2<I2.width() && u2>=0 && 
                v1<I1.height() && v1>=0 && 
                v2<I2.height() && v2>=0  )
            {
               // cout<<"cc1"<<endl;
                c += (I1(u1+x,v1+y) - I1M(u1,v1)) * (I2(u2+x,v2+y) - I2M(u2,v2));
               // cout<<"cc2"<<endl;
        }}}
    return c / ((2*n+1)*(2*n+1));
}

// Compute ZNCC between two patches in images 1 and 2
double zncc(const doubleImage& I1,  // Image 1
            const doubleImage& I1M, // Image of mean intensity value over patch
            const doubleImage& I2,  // Image2
            const doubleImage& I2M, // Image of mean intensity value over patch
            int u1, int v1,         // Pixel of interest in image 1
            int u2, int v2,         // Pixel of interest in image 2
            int n) {                // Half patch size
    double var1 = correl(I1, I1M, I1, I1M, u1, v1, u1, v1, n);
    if (var1 == 0)
        return 0;
    double var2 = correl(I2, I2M, I2, I2M, u2, v2, u2, v2, n);
    if (var2 == 0)
        return 0;
    return correl(I1, I1M, I2, I2M, u1, v1, u2, v2, n) / sqrt(var1 * var2);
}


double p(double& v)
{
    if(v<0){return 1;}
    return sqrt(1-v);
}
// Load two rectified images.
// Compute the disparity of image 2 w.r.t. image 1.
// Display disparity map.
// Display 3D mesh of corresponding depth map.
//
// Images are clipped to focus on pixels visible in both images.
// OPTIMIZATION: to make the program faster, a zoom factor is used to
// down-sample the input images on the fly. For an image of size WxH, you will
// only look at pixels (n+zoom*i,n+zoom*j) with n the radius of patch.
int main() {
    cout << "Loading images... " << flush;
    byteImage I;
    doubleImage I1,I2;
    load(I, srcPath("face00R.png"));
    I1 = I.getSubImage(IntPoint2(20,30), IntPoint2(430,420));
    load(I, srcPath("face01R.png"));
    I2 = I.getSubImage(IntPoint2(20,30), IntPoint2(480,420));
    cout << "done" << endl;

    cout << "Setting parameters... " << flush;
    // Generic parameters
    const int zoom = 2;      // Zoom factor (to speedup computations)
    const int n = 3;         // Consider correlation patches of size (2n+1)*(2n+1)
    const float lambdaf = 0.1; // Weight of regularization (smoothing) term
    const int wcc = max(1+int(1/lambdaf),20); // Energy discretization precision [as we build a graph with 'int' weights]
    const int lambda = lambdaf * wcc; // Weight of regularization (smoothing) term [must be >= 1]
    const float sigma = 3;   // Gaussian blur parameter for disparity
    // Image-specific, hard-coded parameters for approximate 3D reconstruction,
    // as real geometry before rectification is not known
    const int dmin = 10;     // Minmum disparity
    const int dmax = 55;     // Maximum disparity
    const float fB = 40000;  // Depth factor 
    const float db = 100;    // Disparity base
    cout << "done" << endl;

    cout << "Displaying images... " << flush;
    int w1 = I1.width(), w2 = I2.width(), h = I1.height();
    openWindow(w1+w2, h);
    display(grey(I1));
    display(grey(I2), w1, 0);
    cout << "done" << endl;

    cout << "Constructing graph (be patient)... " << flush;
    // Precompute images of mean intensity value over patch
    doubleImage I1M = meanImage(I1,n), I2M = meanImage(I2,n);
    // Zoomed image dimension, disregarding borders (strips of width equal to patch half-size)
    const int nx = (w1 - 2*n) / zoom, ny = (h - 2*n) / zoom;
    const int nd = dmax-dmin; // Disparity range
    const int INF=1000000; // "Infinite" value for edge impossible to cut
    // Create graph
    // The graph library works with node numbers. To clarify the setting, create
    // a formula to associate a unique node number to a triplet (x,y,d) of pixel
    // coordinates and disparity.
    // The library assumes an edge consists of a pair of oriented edges, one in
    // each direction. Put correct weights to the edges, such as 0, INF, or a
    // an intermediate weight.
    /////------------------------------------------------------------
    /////  BEGIN CODE TO COMPLETE: define appropriate graph G
    /////

    Graph<int,int,int> G(nd*nx*ny,0);
    double pzncc = 0, pzncc_last =0;
    double ezncc=0, ezncc_last=0;
    int count = 1;


    int ***nodes = new int**[nd+1];
    for(int i =0; i<nd+1; i++)
    {
        nodes[i] = new int*[w1+1];
        for(int j =0; j<w1+1; j++)
        {
            nodes[i][j] = new int[h+1];
            for(int k = 0; k<h+1;k++)
            {
                nodes[i][j][k] = 0;
            }
        }
    }
    int dminv = -1*dmin;
    
    // int nodes[nd][w1][h];
    cout<<" Constructing graph..."<<endl;
    for(int x=0;x<w1;x++)
    {
        //cout<<"it : ("<<x<<"/"<<w1<<")"<<endl;
        for(int y=0;y<h;y++)
        {
            for(int d=dmin; d<dmax+1; d++)
            {
          //      cout<<"d = "<<d<<endl;
                if(d==dmin)  // FIRST LAYER
                {
                    //cout<<"1"<<endl;
                    nodes[dminv+d][x][y] = G.add_node(1);
                   // cout<<"THis : "<< nodes[dminv+d][x][y]<<endl;
                    //cout<<"2"<<endl;
                    ezncc  = zncc(I1, I1M, I2, I2M, x, y, x+d, y, n);
                    //cout<<"3"<<endl;
                    pzncc = wcc*p(ezncc);
                    //cout<<"4"<<endl;
                    G.add_tweights(nodes[dminv+d][x][y], pzncc, 0);
                    //cout<<"5"<<endl;
                }
                else if(d!=dmax)  // MID LAYERS
                {
                    //cout<<"6"<<endl;
                   //cout<<"non"<<endl;
                   //cout<<"7"<<endl;
                    ezncc  = zncc(I1, I1M, I2, I2M, x, y, x+d, y, n);
                    //cout<<"8"<<endl;
                    pzncc =wcc* p(ezncc);
                    //cout<<"9"<<endl;
                    nodes[dminv+d][x][y] = G.add_node(1); 
                    //cout<<"10"<<endl;
                    G.add_edge(nodes[dminv+d-1][x][y], nodes[dminv+d][x][y], pzncc, 0);
                    //cout<<"11"<<endl;
                }
                else  // LAST LAYER
                {
                    //cout<<"12"<<endl;
                    nodes[dminv+d][x][y] = G.add_node(1); 
                   // cout<<"13"<<endl;
                    ezncc  = zncc(I1, I1M, I2, I2M, x, y, x+d, y, n);
                   // cout<<"14"<<endl;
                    pzncc = wcc*p(ezncc);
                   // cout<<"15"<<endl;
                   //G.add_edge(nodes[dminv+d-1][x][y], nodes[dminv+d][x][y], pzncc, 0);
                   //cout<<"First :"<<nodes[dminv+d-1][x][y]<<endl;
                   //cout<<"Second :"<<nodes[dminv+d][x][y]<<endl;
                   G.add_edge(nodes[dminv+d-1][x][y], nodes[dminv+d][x][y], pzncc, 0);
                   //G.add_edge(29, 30, 28, 0);
                    //cout<<"16"<<endl;
                    ezncc_last  = zncc(I1, I1M, I2, I2M, x, y, x+d+1, y, n);
                    //cout<<"17"<<endl;
                    pzncc_last = wcc*p(ezncc_last);
                    //cout<<"18"<<endl;
                    G.add_tweights(nodes[dminv+d][x][y], 0, pzncc_last);
                    //cout<<"19"<<endl;
                }
                count++;
            }
        }
    }

    cout<<"Adding regularization"<<endl;
    double v = 0;
    for(int x=0;x<w1;x++) // ADDING NEIGHBORS REGULARIZATION
    {
        for(int y=0;y<h;y++)
        {
            for(int d=dmin; d<dmax+1; d++)
            {
                v = lambda;

                if(x-1>=0) G.add_edge(nodes[dminv+d][x-1][y], nodes[dminv+d][x][y], v, v);
                if(x+1<=w1) G.add_edge(nodes[dminv+d][x+1][y], nodes[dminv+d][x][y], v, v);
                if(y-1>=0) G.add_edge(nodes[dminv+d][x][y-1], nodes[dminv+d][x][y], v, v);
                if(y+1<=h) G.add_edge(nodes[dminv+d][x][y+1], nodes[dminv+d][x][y], v, v);
            }
        }
    }

    /* WARNING: dummy code to replace */ 
    /////  END CODE TO BE COMPLETED
    /////------------------------------------------------------------
  
  

  
    cout << "done" << endl;

    cout << "Computing minimum cut... " << flush;
    int f = G.maxflow();
    cout << "done" << endl << "  max flow = " << f << endl;

    cout << "Extracting disparity map from minimum cut... " << flush;
    //doubleImage D(nx,ny);
    doubleImage D(w1, h);
    
    // For each pixel
    for (int i=0;i<w1;i++) 
    {
        for (int j=0;j<h;j++) 
        {
            int d = dmin;
            /*
            cout<<"i = "<<i<<endl;
            cout<<"j = "<<j<<endl;
            cout<<"dminv = "<<dminv<<endl;
            cout<<"d = "<<d<<endl;
            cout<<"1"<<endl;
            */
            while(G.what_segment(nodes[dminv+d][i][j]) == Graph<int,int,int>::SOURCE)
            {
                d += 1;
            }
           // cout<<"2"<<endl;
            D(i, j) = d;
           // cout<<"3"<<endl;
            ///// Extract disparity from minimum cut
            /////------------------------------------------------------------
            /////  BEGIN CODE TO BE COMPLETED: define disparity map D from graph G and minimum cut
            /////
            /* WARNING: dummy code to replace */ // D(i,j) = 40+2*I1(n+i*zoom,n+j*zoom)/(1+(i>j));
            /* WARNING: dummy code to replace */ if (i+j==0) cout<<endl<<"***\n*** Dummy computation: code has to be completed!\n***   (3D generated by default does not make sense)\n***"<<endl;
            /////  END CODE TO BE COMPLETED
            /////------------------------------------------------------------
        }
    }
    cout << "done" << endl;

    cout << "Displaying disparity map... " << flush;
    display(enlarge(grey(D),zoom), n, n);
    cout << "done" << endl;

    cout << "Click to compute and display blured disparity map... " << flush;
    click();
    D=blur(D,sigma);
    display(enlarge(grey(D),zoom),n,n);
    cout << "done" << endl;

    cout << "Click to compute depth map and 3D mesh renderings... " << flush;
    click();
    setActiveWindow( openWindow3D(512, 512, "3D") );
    ///// Compute 3D points
    Array<FloatPoint3> p(nx*ny);
    Array<Color> pcol(nx*ny);
    for(int i=0; i<nx; i++)
        for(int j=0; j<ny; j++) {
            // Compute depth: magic constants depending on camera pose
            float depth = fB / (db + D(i,j)-dmin);
            p[i+nx*j] = FloatPoint3(float(i), float(j), -depth);
            byte g = byte( I1(n+i*zoom,n+j*zoom) );
            pcol[i+nx*j] = Color(g,g,g);
        }
    ///// Create mesh from 3D points
    Array<Triangle> t(2*(nx-1)*(ny-1));
    Array<Color> tcol(2*(nx-1)*(ny-1));
    for (int i=0; i<nx-1; i++)
        for (int j=0; j<ny-1; j++) {
            // Create triangles with next pixels in line/column
            t[2*(i+j*(nx-1))] = Triangle(i+nx*j, i+1+nx*j, i+nx*(j+1));
            t[2*(i+j*(nx-1))+1] = Triangle(i+1+nx*j, i+1+nx*(j+1), i+nx*(j+1));
            tcol[2*(i+j*(nx-1))] = pcol[i+nx*j];
            tcol[2*(i+j*(nx-1))+1] = pcol[i+nx*j];
        }
    // Create first mesh as textured with colors taken from original image
    Mesh Mt(p.data(), nx*ny, t.data(), 2*(nx-1)*(ny-1), 0, 0, FACE_COLOR);
    Mt.setColors(TRIANGLE, tcol.data());
    // Create second mesh with artificial light
    Mesh Mg(p.data(), nx*ny, t.data(), 2*(nx-1)*(ny-1), 0, 0,
            CONSTANT_COLOR, SMOOTH_SHADING);
    cout << "done" << endl;

    // Display 3D mesh renderings
    cout << "***** 3D mesh renderings *****" << endl;
    cout << "- Button 1: toggle textured or gray rendering" << endl;
    cout << "- SHIFT+Button 1: rotate" << endl;
    cout << "- SHIFT+Button 3: translate" << endl;
    cout << "- Mouse wheel: zoom" << endl;
    cout << "- SHIFT+a: zoom out" << endl;
    cout << "- SHIFT+z: zoom in" << endl;
    cout << "- SHIFT+r: recenter camera" << endl;
    cout << "- SHIFT+m: toggle solid/wire/points mode" << endl;
    cout << "- Button 3: exit" << endl;
    showMesh(Mt);
    bool textured = true;
    while (true) {
        Event evt;
        getEvent(5, evt);
        // On mouse button 1
        if (evt.type == EVT_BUT_ON && evt.button == 1) {
            // Toggle textured rendering and gray rendering
            if (textured) {
                hideMesh(Mt,false);
                showMesh(Mg,false);
            } else {
                hideMesh(Mg,false);
                showMesh(Mt,false);
            }
            textured = !textured;
        }
        // On mouse button 3
        if (evt.type == EVT_BUT_ON && evt.button == 3)
            break;
    }

    endGraphics();
    return 0;
}
