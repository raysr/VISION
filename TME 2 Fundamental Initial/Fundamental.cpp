// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse
// Date:     2013/10/08

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include<cmath>
#include <stdlib.h> 
#include <algorithm>

using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    cout.precision(2000);
    
    
    float beta = 0.01;
    
    
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS

    // SVD
   
   // NORMALIZATION  
    FMatrix<float, 3, 3> N; 
    N(0, 0) = 0.001; N(0, 1) = 0; N(0, 2) = 0;
    N(1, 0) = 0; N(1, 1) = 0.001; N(1, 2) = 0; 
    N(2, 0) = 0; N(2, 1) = 0; N(2, 2) = 1;
    srand (time(NULL));
    int best_performance = 0;
    vector<int> bestInliers;
    float niter = 1000;

    for(int lrm=0; lrm<niter; lrm++) // RANSAC LOOP
    {
        
        cout<<"("<<lrm<<"/"<<niter<<")"<<endl;


    // cout<<"T1"<<endl;
    int size = 8;
    FMatrix<float, 8, 9> A;
    FVector<float, 3> v1, v2;
    vector<Match> matchestmp(8);
    int randnu;
    for(int i=0; i<size; i++)
    {
        randnu =  rand() % matches.size() + 0;
        v1[0] = matches.at(randnu).x1; v1[1] = matches.at(randnu).y1 ; v1[2] = 1;
        v2[0] = matches.at(randnu).x2; v2[1] = matches.at(randnu).y2 ; v2[2] = 1;
        v1 = N*v1;
        v2 = N*v2;

        matchestmp[i].x1 = v1[0];
        matchestmp[i].y1 = v1[1];
        matchestmp[i].x2 = v2[0];
        matchestmp[i].y2 = v2[1];
        // cout<<"before : "<<matches.at(i).x1<<","<<matches.at(i).y1<<", "<<matches.at(i).x2<<", "<<matches.at(i).y2<<endl;
        // cout<<"after : "<<matchestmp[i]<<endl;
    }


     //cout<<"T2"<<endl;
    for(int i=0; i<size; i++)
    {
        A(i, 0) = matchestmp[i].x2*matchestmp[i].x1;
        A(i, 1) = matchestmp[i].x2*matchestmp[i].y1;
        A(i, 2) = matchestmp[i].x2;
        A(i, 3) = matchestmp[i].y2*matchestmp[i].x1;
        A(i, 4) = matchestmp[i].y2*matchestmp[i].y1;
        A(i, 5) = matchestmp[i].y2;
        A(i, 6) = matchestmp[i].x1;
        A(i, 7) = matchestmp[i].y1;
        A(i, 8) = 1;

    }

    // cout<<"Matrix A :"<<endl;
    // cout<<A<<endl;

    // A = U * D * V
    FMatrix<float, 8, 8> U;
    FVector<float, 8> D;
    FMatrix<float, 9, 9> Vt;
    svd(A, U, D, Vt) ;  
    /*
    cout<<"U ="<<endl;
    cout<<U;
    cout<<endl;
    cout<<"D ="<<endl;
    cout<<D;
    cout<<endl;
    cout<<"V^T ="<<endl;
    cout<<Vt;
    cout<<endl;
    */

    FVector<float, 9> f = Vt.getRow(8);
    
    /*
     cout<<"f = V(8,:) ="<<endl;
     cout<<f<<endl;
    */


    // reshaping f into F 
    // cout<<"T3"<<endl;
    FMatrix<float, 3, 3> F;
    F(0, 0) = f[0]; F(0, 1) = f[1]; F(0, 2) = f[2];
    F(1, 0) = f[3]; F(1, 1) = f[4]; F(1, 2) = f[5];
    F(2, 0) = f[6]; F(2, 1) = f[7]; F(2, 2) = f[8]; 


   // cout<<"RESULTING F :"<<endl;
   // cout<<F<<endl;

    //  Reducing Rank to 2 by setting lowest value of Df to 0
    FMatrix<float, 3, 3> Uf;
    FVector<float, 3> Df;
    FMatrix<float, 3, 3> Vft;
    svd(F, Uf, Df, Vft) ;  
     /*
     cout<<"Uf ="<<endl;
     cout<<Uf;
     cout<<endl;
     cout<<"Df ="<<endl;
     cout<<Df;
     cout<<endl;
     cout<<"Vf^T ="<<endl;
     cout<<Vft;
     cout<<endl;
    */

    Df[2] = 0;
    FMatrix<float, 3, 3> tmp = Diagonal(Df);
    F = Uf*tmp*Vft;

    // Last step of Normalization
    F = transpose(N)*F*N;
    float distance;
    long double number_inliers=0, number_outliers=0;
    float treshold = 0.1;
   //  cout<<"T4"<<endl;
    int pr = matches.size();
    int training_size = min(pr, 100);


    vector<int> Inliers;
    for(int i=0;i<training_size;i++)
    {
        // SAMPSON DISTANCE
        //  cout<<"i :"<<i<<endl;
        // cout<<"T4.1"<<endl;
        
        v1[0] = matches.at(i).x1; v1[1] = matches.at(i).y1 ; v1[2] = 1;
        v2[0] = matches.at(i).x2; v2[1] = matches.at(i).y2 ; v2[2] = 1;
        // cout<<"T4.2"<<endl;
        FVector<float, 3> tmp = transpose(F)*v1;
        // cout<<"T4.3"<<endl;
        
        distance = abs( tmp[0]*v2[0]   +   tmp[1]*v2[1]  +    tmp[2] )   /   sqrt( pow(tmp[0], 2.0 )  +   pow( tmp[1], 2.0 )  );
        // cout<<"T4.4"<<endl;
        // cout<<"Distance = "<<distance<<endl;
        if(distance<treshold)
        {
            // cout<<"T4.4.A1"<<endl;
            Inliers.push_back(i);
            // cout<<"T4.4.A2"<<endl;
            number_inliers++;
            // cout<<"T4.4.A3"<<endl;
        }
        else
        {
            // cout<<"T4.4.B1"<<endl;
            number_outliers++;
            // cout<<"T4.4.B2"<<endl;
        }
        // cout<<"T4.5"<<endl;
    }
    // cout<<"T5"<<endl;
    
    cout<<"Number of inliers = "<<number_inliers<<endl;
    cout<<"Number of outliers = "<<number_outliers<<endl;
    
    if(number_inliers>best_performance)
    {
        best_performance = number_inliers;
        long double a = number_inliers/(number_inliers+number_outliers);
        long double b = pow((a), 8.0);
        long double tmp = 1.0 - b;
        long double c = log(tmp);
       /*
        cout<<"log(beta) = "<<log(beta)<<endl;
        cout<<"number_inliers/(number_inliers+number_outliers) ="<<a<<endl;
        cout<<"(num inliers/sum)^8 = "<<b<<endl;
        cout<<"1-b = "<<tmp<<endl;
        cout<<"log(1-pow((number_inliers/(number_inliers+number_outliers)), 8)) = "<<c<<endl;
        */
        niter = log(beta)/c;
        FMatrix<float,3,3> bestF(F);
        vector<int>  bestInliers(Inliers);
    }
    //cout<<"T6"<<endl;
    
    }
    cout<<"BEST PERFORMANCE NUMBER INLIERS = "<<best_performance<<endl;
    
    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;
        // --------------- TODO ------------
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    cout << " matches: " << matches.size() << endl;
    click();
    
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
