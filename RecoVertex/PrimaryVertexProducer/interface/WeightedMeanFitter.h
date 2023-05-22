#ifndef RecoVertex_PrimaryVertexProducer_WeightedMeanFitter_h
#define RecoVertex_PrimaryVertexProducer_WeightedMeanFitter_h

#include <vector>
#include <math.h> // JS_EDIT
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#define DEBUG
#ifdef DEBUG
#define DEBUGLEVEL -5 // -5 means: debug the same way as CPU; Final value TBD
#endif

namespace WeightedMeanFitter {

  constexpr double startError = 20.0;
  constexpr double precision = 1e-24;
  constexpr double corr_x = 1.2;
  constexpr double corr_x_bs = 1.0; // corr_x for beam spot 
  constexpr double corr_z = 1.4;
  constexpr int maxIterations = 50;

  //inline std::pair<GlobalPoint, double> nearestPoint(const GlobalPoint& vertex, reco::Track iclus){
  inline std::pair<GlobalPointDouble, double> nearestPoint(const GlobalPointDouble& vertex, reco::TransientTrack iclus){

#ifdef DEBUG
    printf("cpu fitter vx,px,vertexx %.10f,%.10f,%.10f \n",iclus.stateAtBeamLine().trackStateAtPCA().position().x(),iclus.track().px(),vertex.x());
#endif

      double ox = iclus.stateAtBeamLine().trackStateAtPCA().position().x();
      double oy = iclus.stateAtBeamLine().trackStateAtPCA().position().y();
      double oz = iclus.stateAtBeamLine().trackStateAtPCA().position().z();

      double vx = iclus.track().px();
      double vy = iclus.track().py();
      double vz = iclus.track().pz();

      double opx = vertex.x() - ox;
      double opy = vertex.y() - oy;
      double opz = vertex.z() - oz;

      double vnorm2 = (vx*vx + vy*vy + vz*vz);
      double t = (vx * opx + vy * opy + vz * opz) / (vnorm2);

      GlobalPointDouble p(ox + t * vx, oy + t * vy, oz + t * vz);
      return std::pair<GlobalPointDouble, double>(p, std::sqrt( std::pow( p.x() - vertex.x() , 2) + std::pow( p.y() - vertex.y() , 2) + std::pow( p.z() - vertex.z() , 2) ));
  }

  inline TransientVertex weightedMeanOutlierRejection(const std::vector<std::pair<GlobalPointDouble, GlobalPointDouble>>& points, std::vector<reco::TransientTrack> iclus){
       double x=0., y=0., z=0.;
       double s_wx=0., s_wz=0.;
       double s2_wx=0., s2_wz=0.;
       double wx=0., wz=0., chi2=0.;
       double ndof_x = 0.;

       AlgebraicSymMatrix33 err;
       err(0,0) = startError/10 * startError/10;
       err(1,1) = startError/10 * startError/10;
       err(2,2) = startError * startError; // error is 20 cm, so cov -> is 20 ^ 2
       for (const auto& p : points){

              wx = p.second.x() <=  precision ? 1. / std::pow(precision,2) : 1. / std::pow(p.second.x(),2);

              wz = p.second.z() <=  precision ? 1. / std::pow(precision,2) : 1. / std::pow(p.second.z(),2);

             x += p.first.x() * wx;
             y += p.first.y() * wx;
             z += p.first.z() * wz;

             s_wx += wx;
             s_wz += wz;

#ifdef DEBUG
	     printf("cpu fitter x,dx2,wx,z,dz2,wz %.10f,%.10f,%.10f,%.10f,%.10f,%.10f \n",p.first.x(),p.second.x()*p.second.x(),wx,p.first.z(),p.second.z()*p.second.z(),wz);
#endif
       }

       if ( s_wx == 0. || s_wz == 0. ){
          edm::LogWarning("WeightedMeanFitter") << "Vertex fitting failed at beginning \n";
          return TransientVertex(GlobalPointDouble(0,0,0), err, iclus, 0, 0);
       }

       x /= s_wx;
       y /= s_wx;
       z /= s_wz;

      double old_x, old_y, old_z;

      double xpull;
      int niter = 0;
      double mu = 3.;

      double err_x, err_z;

      err_x = 1. / s_wx;
      err_z = 1. / s_wz;

#ifdef DEBUG
    printf("cpu fitter x,err_x,z,err_z %.10f,%.10f,%.10f,%.10f \n",x,err_x,z,err_z);
#endif

      while ((niter++) < 2){
          old_x = x;
          old_y = y;
          old_z = z;
          s_wx = 0; s_wz = 0; s2_wx = 0; s2_wz = 0;
          x = 0; y = 0; z = 0;
          ndof_x = 0;

          for (unsigned int i = 0; i < (unsigned int) points.size(); i++){

	      //std::pair<GlobalPoint, double> p = nearestPoint(GlobalPoint(old_x, old_y, old_z), (iclus)[i].track());
              std::pair<GlobalPointDouble, double> p = nearestPoint(GlobalPointDouble(old_x, old_y, old_z), (iclus)[i]);

              wx =  points[i].second.x() <= precision ? std::pow(precision, 2) : std::pow(points[i].second.x(), 2);
              wz =  points[i].second.z() <= precision ? std::pow(precision, 2) : std::pow(points[i].second.z(),2);

#ifdef DEBUG
	      printf("cpu fitter niter,x,dx2,wx,z,dz2,wz %d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f \n",niter,p.first.x(),points[i].second.x()*points[i].second.x(),wx,p.first.z(),points[i].second.z()*points[i].second.z(),wz);
#endif
              xpull = 0.;

              if ( std::pow(p.first.x() - old_x, 2) / (wx + err_x) < mu*mu
                && std::pow(p.first.y() - old_y, 2) / (wx + err_x) < mu*mu
                && std::pow(p.first.z() - old_z, 2) / (wz + err_z) < mu*mu) {
                  xpull = 1.;
#ifdef DEBUG
		  printf("cpu fitter: adding track");
#endif

              }
              //// gaussian weighting
	      /*
              if (xpull == 1) {
                  double coeff = std::pow(2/(M_PI*(wz+err_z)),0.25);
                  double func = std::exp(-1*std::pow(p.first.z() - old_z, 2) / (wz + err_z));
                  xpull = coeff * func; //creates the idea of "influence". take off coeff for non-normalized
              }
	      */

              ndof_x += xpull;

              wx = xpull / wx;
              wz = xpull / wz;

              x += wx * p.first.x();
              y += wx * p.first.y();
              z += wz * p.first.z();

              s_wx += wx;
              s_wz += wz;

              s2_wx += wx * xpull;
              s2_wz += wz * xpull;
          }

          if ( s_wx == 0. || s_wz == 0. ){
              edm::LogWarning("WeightedMeanFitter") << "Vertex fitting failed" << s_wx << " , " << s_wz << "\n";
              return TransientVertex(GlobalPointDouble(0,0,0), err, iclus, 0, 0);
          }
          x /= s_wx;
          y /= s_wx;
          z /= s_wz;

          err_x = (s2_wx / std::pow(s_wx, 2));
          err_z = (s2_wz / std::pow(s_wz, 2));

          if (std::abs(x - old_x) < (precision/1.) && std::abs(y - old_y) < (precision/1.) && std::abs(z - old_z) < (precision/1.)){
              break;
          }
       } 
       
       err(0,0) = err_x * corr_x * corr_x;
       err(1,1) = err_x * corr_x * corr_x;
       err(2,2) = err_z * corr_z * corr_z;


       double dist = 0;
       for (const auto& p : points){
          wx = p.second.x();
          wx =  wx <= precision ? precision : wx;

          wz = p.second.z();
          wz =  wz <= precision ? precision : wz;

          dist =  std::pow(p.first.x() - x, 2) / ( std::pow(wx, 2) +  err(0,0) );
          dist += std::pow(p.first.y() - y, 2) / ( std::pow(wx, 2) +  err(1,1) );
          dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) +  err(2,2) );
          chi2 += dist;
       }

#ifdef DEBUG
       printf("end of cpu fitter x,y,z %.10f,%.10f,%.10f \n",x,y,z);
#endif

       TransientVertex v(GlobalPointDouble(x,y,z), err, iclus, chi2, (int) ndof_x);
       return v;
  }



inline TransientVertex weightedMeanOutlierRejectionBeamSpot(const std::vector<std::pair<GlobalPointDouble, GlobalPointDouble>>& points, std::vector<reco::TransientTrack> iclus, const reco::BeamSpot& beamSpot){
     double x=0., y=0., z=0.;
     double s_wx=0., s_wz=0.;
     double s2_wx=0., s2_wz=0.;
     double wx=0., wz=0., chi2=0.;
     double wy=0., s_wy=0., s2_wy=0.;
     double ndof_x = 0.; 

     AlgebraicSymMatrix33 err;
     err(0,0) = startError/10 * startError/10;
     err(1,1) = startError/10 * startError/10;
     err(2,2) = startError * startError; // error is 20 cm, so cov -> is 20 ^ 2

     GlobalError bse(beamSpot.rotatedCovariance3D());
     GlobalPointDouble bsp(Basic3DVector<double>(beamSpot.position()));
   
     for (const auto& p : points){ 

            wx = p.second.x() <=  precision ? 1. / std::pow(precision,2) : 1. / std::pow(p.second.x(),2);
            wy = p.second.y() <=  precision ? 1. / std::pow(precision,2) : 1. / std::pow(p.second.y(),2);

            wz = p.second.z() <=  precision ? 1. / std::pow(precision,2) : 1. / std::pow(p.second.z(),2);

           x += p.first.x() * wx;
           y += p.first.y() * wy;
           z += p.first.z() * wz;

           s_wx += wx;
           s_wy += wy;
           s_wz += wz;
	   
#ifdef DEBUG
	   printf("cpu fitter x,dx2,wx,z,dz2,wz %.10f,%.10f,%.10f,%.10f,%.10f,%.10f \n",p.first.x(),p.second.x()*p.second.x(),wx,p.first.z(),p.second.z()*p.second.z(),wz);
#endif

     }

     if ( s_wx == 0. || s_wy == 0. ||s_wz == 0. ){
        edm::LogWarning("WeightedMeanFitter") << "Vertex fitting failed at beginning \n";
        return TransientVertex(GlobalPointDouble(0,0,0), err, iclus, 0, 0);
     }

#ifdef DEBUG
        printf("cpu fitter using beamspot x,cxx,y,cyy %.10f,%.10f,%.10f,%.10f \n",bsp.x(),bse.cxx(),bsp.y(),bse.cyy());
#endif

    // use the square of covariance element to increase it's weight: it will be the most important 
     wx = bse.cxx() <=  precision ? 1. / std::pow(precision,2) : 1. / std::pow(bse.cxx(),2);
     wy = bse.cyy() <=  precision ? 1. / std::pow(precision,2) : 1. / std::pow(bse.cyy(),2);

     x += bsp.x() * wx;
     y += bsp.y() * wy;
     
     x /= (s_wx + wx);     
     y /= (s_wy + wy);     
     z /= s_wz;  
 
    double old_x, old_y, old_z;

    double xpull;
    int niter = 0;
    double mu = 3.;

    double err_x, err_y, err_z;

    err_x = 1. / s_wx;
    err_y = 1. / s_wy;
    err_z = 1. / s_wz;

#ifdef DEBUG
    printf("cpu fitter x,err_x,z,err_z %.10f,%.10f,%.10f,%.10f \n",x,err_x,z,err_z);
#endif

    while ((niter++) < 2){
        old_x = x;
        old_y = y;
        old_z = z;
        s_wx = 0; s_wz = 0; s2_wx = 0; s2_wz = 0;

        s_wy = 0; s2_wy = 0;

        x = 0; y = 0; z = 0;
        ndof_x = 0; 

        for (unsigned int i = 0; i < (unsigned int) points.size(); i++){ 
	    //std::pair<GlobalPoint, double> p = nearestPoint(GlobalPoint(old_x, old_y, old_z), (iclus)[i].track());
            std::pair<GlobalPointDouble, double> p = nearestPoint(GlobalPointDouble(old_x, old_y, old_z), (iclus)[i]);


            wx =  points[i].second.x() <= precision ? std::pow(precision, 2) : std::pow(points[i].second.x(), 2);
            wy =  points[i].second.y() <= precision ? std::pow(precision, 2) : std::pow(points[i].second.y(), 2);

            wz =  points[i].second.z() <= precision ? std::pow(precision, 2) : std::pow(points[i].second.z(),2);
           
#ifdef DEBUG
	      printf("cpu fitter niter,x,dx2,wx,z,dz2,wz %d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f \n",niter,p.first.x(),points[i].second.x()*points[i].second.x(),wx,p.first.z(),points[i].second.z()*points[i].second.z(),wz);
#endif
 
            xpull = 0.;
            if ( std::pow(p.first.x() - old_x, 2) / (wx + err_x) < mu*mu && std::pow(p.first.y() - old_y, 2) / (wy + err_y) < mu*mu && std::pow(p.first.z() - old_z, 2) / (wz + err_z) < mu*mu)  xpull = 1.;

            ndof_x += xpull;

            wx = xpull / wx;
            wy = xpull / wy;
            wz = xpull / wz;

            x += wx * p.first.x();
            y += wy * p.first.y();
            z += wz * p.first.z();

            s_wx += wx;
            s_wy += wy;
            s_wz += wz;

            s2_wx += wx * xpull;
            s2_wy += wy * xpull;
            s2_wz += wz * xpull;
        }

        if ( s_wx == 0. || s_wy == 0. || s_wz == 0. ){
            edm::LogWarning("WeightedMeanFitter") << "Vertex fitting failed" << s_wx << " , " << s_wy << " , " << s_wz << "\n";
            return TransientVertex(GlobalPointDouble(0,0,0), err, iclus, 0, 0);
        }

        wx = bse.cxx() <=  std::pow(precision,2) ? 1. / std::pow(precision,2) : 1. / bse.cxx();
        wy = bse.cyy() <=  std::pow(precision,2) ? 1. / std::pow(precision,2) : 1. / bse.cyy();

        x += bsp.x() * wx;
        y += bsp.y() * wy;
        s_wx  += wx;
        s2_wx += wx;
        s_wy  += wy;
        s2_wy += wy;
        
        x /= s_wx;     
        y /= s_wy;     
        z /= s_wz;

        err_x = (s2_wx / std::pow(s_wx, 2));
        err_y = (s2_wy / std::pow(s_wy, 2));
        err_z = (s2_wz / std::pow(s_wz, 2));

        if (std::abs(x - old_x) < (precision/1.) && std::abs(y - old_y) < (precision/1.) && std::abs(z - old_z) < (precision/1.)){
            break;
        }
     }
     err(0,0) = err_x * corr_x_bs * corr_x_bs;
     err(1,1) = err_y * corr_x_bs * corr_x_bs;
     err(2,2) = err_z * corr_z * corr_z;


     double dist = 0; 
     for (const auto& p : points){ 
        wx = p.second.x();
        wx =  wx <= precision ? precision : wx;

        wz = p.second.z();
        wz =  wz <= precision ? precision : wz;

        dist =  std::pow(p.first.x() - x, 2) / ( std::pow(wx, 2) +  err(0,0) );
        dist += std::pow(p.first.y() - y, 2) / ( std::pow(wx, 2) +  err(1,1) );
        dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) +  err(2,2) ); 
        chi2 += dist;

     }

#ifdef DEBUG
     printf("end of cpu fitter x,y,z %.10f,%.10f,%.10f \n",x,y,z);
#endif

     TransientVertex v(GlobalPointDouble(x,y,z), err, iclus, chi2, (int) ndof_x);
     return v;
}

  inline TransientVertex weightedMeanOutlierRejectionVarianceAsError(const std::vector<std::pair<GlobalPointDouble, GlobalPointDouble>>& points, std::vector<std::vector<reco::TransientTrack>>::const_iterator iclus){
       double x=0, y=0, z=0, s_wx=0, s_wy=0, s_wz=0, s2_wx=0, s2_wy=0, s2_wz=0, wx=0, wy=0, wz=0, chi2=0;
       double ndof_x = 0;
       AlgebraicSymMatrix33 err;
       err(0,0) = startError/10 * startError/10;
       err(1,1) = startError/10 * startError/10;
       err(2,2) = startError * startError; // error is 20 cm, so cov -> is 20 ^ 2

       for (const auto& p : points){
              wx = p.second.x();
              wx = wx <=  precision ? 1. / std::pow(precision,2) : 1. / std::pow(wx,2);

              wz = p.second.z();
              wz = wz <=  precision ? 1. / std::pow(precision,2) : 1. / std::pow(wz,2);

              x += p.first.x() * wx;
              y += p.first.y() * wx;
              z += p.first.z() * wz;

              s_wx += wx;
              s_wz += wz;
       }

       if ( s_wx == 0. || s_wz == 0. ){
          edm::LogWarning("WeightedMeanFitter") << "Vertex fitting failed at beginning. \n";
          return TransientVertex(GlobalPointDouble(0,0,0), err, *iclus, 0, 0);
       }

       x /= s_wx;
       y /= s_wx;
       z /= s_wz;

      double old_x, old_y, old_z;
      double xpull;
      int niter = 0;
      double mu = 3.;
      double err_x, err_y, err_z;
      err_x = 1. / std::sqrt(s_wx);
      err_y = 1. / std::sqrt(s_wx);
      err_z = 1. / std::sqrt(s_wz);
      double s_err_x = 0, s_err_y = 0, s_err_z = 0;
      while ((niter++) < maxIterations){
          old_x = x;
          old_y = y;
          old_z = z;
          s_wx = 0; s_wy = 0; s_wz = 0; s2_wx = 0; s2_wy = 0; s2_wz = 0;
          x = 0; y = 0; z = 0;
          s_err_x = 0.; s_err_y = 0.; s_err_z = 0;

          for (const auto& p : points){
              wx = p.second.x();
              wx =  wx <= precision ? precision : wx;

              wy = wx*wx + err_y*err_y;
              wx = wx*wx + err_x*err_x;

              wz = p.second.z();
              wz =  wz <= precision ? precision : wz;
              wz = wz*wz + err_z*err_z;

              xpull  = std::pow((p.first.x() - old_x), 2) / wx;
              xpull += std::pow((p.first.y() - old_y), 2) / wy;
              xpull += std::pow((p.first.z() - old_z), 2) / wz;
              xpull = 1. / (1. + std::exp(-0.5 * ((mu * mu) - xpull)));
              ndof_x += xpull;

              wx = 1. / wx;
              wy = 1. / wy;
              wz = 1. / wz;

              wx *= xpull;
              wy *= xpull;
              wz *= xpull;

              x += wx * p.first.x();
              y += wy * p.first.y();
              z += wz * p.first.z();

              s_wx += wx;
              s_wy += wy;
              s_wz += wz;

              s2_wx += wx * wx;
              s2_wy += wy * wy;
              s2_wz += wz * wz;

            s_err_x += wx * pow(p.first.x() - old_x, 2);
            s_err_y += wy * pow(p.first.y() - old_y, 2);
            s_err_z += wz * pow(p.first.z() - old_z, 2);
          }
          if ( s_wx == 0. || s_wy == 0. || s_wz == 0. ){

            edm::LogWarning("WeightedMeanFitter") << "Vertex fitting failed" << s_wx << " , " << s_wy << " , " << s_wz << "\n";
              return TransientVertex(GlobalPointDouble(0,0,0), err, *iclus, 0, 0);
          }
          x /= s_wx;
          y /= s_wy;
          z /= s_wz;

          err_x = std::sqrt(s_err_x / s_wx );
          err_y = std::sqrt(s_err_y / s_wy );
          err_z = std::sqrt(s_err_z / s_wz );

          if (std::abs(x - old_x) < (precision/1.) && std::abs(y - old_y) < (precision/1.) && std::abs(z - old_z) < (precision/1.))
              break;
      }

       err(0,0) = err_x * err_x;
       err(1,1) = err_y * err_y;
       err(2,2) = err_z * err_z;

       double dist = 0;
       for (const auto& p : points){
          wx = p.second.x();
          wx =  wx <= precision ? precision : wx;

          wz = p.second.z();
          wz =  wz <= precision ? precision : wz;

          dist =  std::pow(p.first.x() - x, 2) / ( std::pow(wx, 2) + std::pow( err(0,0), 2) );
          dist += std::pow(p.first.y() - y, 2) / ( std::pow(wx, 2) + std::pow( err(1,1), 2) );
          dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) + std::pow( err(2,2), 2) );
          chi2 += dist;
       }
       TransientVertex v(GlobalPointDouble(x,y,z), err, *iclus, chi2, (int) ndof_x);
       return v;
  }

};

#endif
