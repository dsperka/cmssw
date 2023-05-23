// CUDA include files
#include <cuda_runtime.h>

// CMSSW include files
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "RecoVertex/PrimaryVertexProducer/interface/fitterCUDA.h"
#include "CUDADataFormats/Track/interface/TrackForPVHeterogeneous.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
//#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
//#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerCUDA.h" //JS
//#include "RecoVertex/PrimaryVertexProducer/interface/WeightedMeanFitter.h" //JS
//#include "DataFormats/Math/interface/Error.h"
//#include <cstddef>
//#include <cstdint>
#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"
//#include <thrust/device_vector.h>
//#include <thrust/host_vector.h>
//#include <thrust/sort.h>
//#include "HeterogeneousCore/CUDAUtilities/interface/radixSort.h"
#include <math.h>

//#define DEBUG

using Vector512d = Eigen::Matrix<double, 1024, 1>;

namespace fitterCUDA {

__global__ void fitterKernel(
    unsigned int ntracks,
    TrackForPV::TrackForPVSoA* tracks,
    TrackForPV::VertexForPVSoA* vertices,
    algo algorithm,
    bs beamspot
){

  size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x;
  size_t gridSize = blockDim.x * gridDim.x;
  double precision = 1e-24; //FLOAT CHANGE
  double precision2 = precision*precision;

  //1 block for each vertex (for block in BLOCKS)
  for (unsigned int k = firstElement; k < vertices->nTrueVertex(0); k += gridSize) {
    unsigned int ivertex = vertices->order(k);
    if (!vertices->isGood(ivertex)) continue; //skip if not good

#ifdef DEBUG
    printf("start gpu fitter k=%d,ivtx=%d,ntracks=%d,z=%.10f \n",k,ivertex,vertices->ntracks(ivertex),vertices->z(ivertex));
#endif

    //-----------------------------position/vector loop-----------------------------
    //1 thread / track (for each thread in block)
    //double err_x = 4.0, err_y = 4.0, err_z = 400.0; //this looks like (startError/10)^2 in CPU code
    double corr_x = 1.2, corr_x_bs = 1.0, corr_y = 1.2, corr_z = 1.4;
    Vector512d track_ids, track_weights;
    int track_id_counter = 0, track_weight_counter = 0;

    //loop to calculate z weighted average
    double x=0., y=0., z=0.;
    double s_wx=0., s_wz=0.;
    double s2_wx=0., s2_wz=0.;
    double wx=0., wz=0., chi2=0.;
    double wy=0., s_wy=0., s2_wy=0.;
    vertices->ndof(ivertex) = 0; 

    for (unsigned int kk = 0; kk < tracks->nTrueTracks; kk++){
      unsigned int itrack = tracks->order(kk);
      if (not(tracks->isGood(itrack))) continue;
      unsigned int ivtxFromTk = tracks->kmin(itrack);

      if (ivtxFromTk == ivertex) {

        track_ids[track_id_counter] = itrack;
        track_id_counter++;

	if (algorithm.useBeamConstraint) {
	  wx = (tracks->dxy2(itrack) <= precision2) ? 1./(precision2) : 1./tracks->dxy2(itrack);
	  wy = (tracks->dxy2(itrack) <= precision2) ? 1./(precision2) : 1./tracks->dxy2(itrack);
	  wz = (tracks->dz2(itrack) <= precision2) ? 1./(precision2) : 1./tracks->dz2(itrack);

	  x += tracks->x(itrack) * wx; 
	  y += tracks->y(itrack) * wy;
	  z += tracks->z(itrack) * wz;

	  s_wx += wx;
	  s_wy += wy;
	  s_wz += wz;

#ifdef DEBUG
	printf("gpu fitter x,dx2,wx,z,dz2,wz %.10f,%.10f,%.10f,%.10f,%.10f,%.10f \n",tracks->x(itrack),tracks->dxy2(itrack),wx,tracks->z(itrack),tracks->dz2(itrack),wz);
#endif

	} else {
	  wx = (tracks->dxy2AtIP(itrack) <= precision2) ? 1./(precision2) : 1./tracks->dxy2AtIP(itrack);
	  wy = (tracks->dxy2AtIP(itrack) <= precision2) ? 1./(precision2) : 1./tracks->dxy2AtIP(itrack);
	  wz = (tracks->dz2(itrack) <= precision2) ? 1./(precision2) : 1./tracks->dz2(itrack);

	  x += tracks->xAtIP(itrack) * wx; 
 	  y += tracks->yAtIP(itrack) * wy;
	  z += tracks->z(itrack) * wz;

	  s_wx += wx;
	  s_wy += wy;
	  s_wz += wz;

#ifdef DEBUG
	printf("gpu fitter x,dx2,wx,z,dz2,wz %.10f,%.10f,%.10f,%.10f,%.10f,%.10f \n",tracks->xAtIP(itrack),tracks->dxy2AtIP(itrack),wx,tracks->z(itrack),tracks->dz2(itrack),wz);
#endif

	}
      }
    }

    vertices->track_id(ivertex) = track_ids;

    if (algorithm.useBeamConstraint) {

#ifdef DEBUG
	printf("gpu fitter using beamspot x,cxx,y,cyy %.10f,%.10f,%.10f,%.10f \n",beamspot.x,beamspot.cxx,beamspot.y,beamspot.cyy);
#endif

      wx = beamspot.cxx <=  precision ? 1./precision2 : 1. / (beamspot.cxx*beamspot.cxx);
      wy = beamspot.cyy <=  precision ? 1./precision2 : 1. / (beamspot.cyy*beamspot.cyy);
    
      x += beamspot.x * wx;
      y += beamspot.y * wy;
      
      x /= (s_wx + wx);
      y /= (s_wy + wy);
      z /= s_wz;

    } else {

      x /= s_wx;
      y /= s_wy;
      z /= s_wz;

    }

    double old_x, old_y, old_z;

    double xpull;
    int niter = 0;
    double mu = 3.;
    double mu2 = mu*mu;

    double err_x, err_y, err_z;

    err_x = 1. / s_wx;
    err_y = 1. / s_wy;
    err_z = 1. / s_wz;

#ifdef DEBUG
    printf("gpu fitter x,err_x,z,err_z %.10f,%.10f,%.10f,%.10f \n",x,err_x,z,err_z);
#endif

    while ((niter++) < 2){
      old_x = x;
      old_y = y;
      old_z = z;

      s_wx = 0; s_wy = 0; s_wz = 0;
      s2_wx = 0; s2_wy = 0; s2_wz = 0;

      x = 0; y = 0; z = 0;
      vertices->ndof(ivertex) = 0;

      for (unsigned int kk = 0; kk < tracks->nTrueTracks; kk++){

	unsigned int itrack = tracks->order(kk);
        if (not(tracks->isGood(itrack))) continue;

	unsigned int ivtxFromTk = tracks->kmin(itrack);
	if (ivtxFromTk == ivertex) {
	  
	  double ox = tracks->x(itrack); 
	  double oy = tracks->y(itrack); 
	  double oz = tracks->z(itrack);

	  double vx = tracks->px(itrack);
	  double vy = tracks->py(itrack);
	  double vz = tracks->pz(itrack);

#ifdef DEBUG
	  printf("gpu fitter vx,px,vertexx %.10f,%.10f,%.10f \n",ox,vx,old_x);
#endif

	  double opx = old_x - ox;
	  double opy = old_y - oy;
	  double opz = old_z - oz;

	  double vnorm2 = (vx*vx + vy*vy + vz*vz);
	  double t = (vx * opx + vy * opy + vz * opz) / (vnorm2);

	  double tx = ox + t * vx;
	  double ty = oy + t * vy;
	  double tz = oz + t * vz;

	  ////WEIGHTING
	  double wx, wy, wz;

	  if (algorithm.useBeamConstraint) {
	    wx = (tracks->dxy2(itrack) <= precision2) ? precision2 : tracks->dxy2(itrack);
	    wy = (tracks->dxy2(itrack) <= precision2) ? precision2 : tracks->dxy2(itrack);
	  } else {
	    wx = (tracks->dxy2AtIP(itrack) <= precision2) ? precision2 : tracks->dxy2AtIP(itrack);
	    wy = (tracks->dxy2AtIP(itrack) <= precision2) ? precision2 : tracks->dxy2AtIP(itrack);
	  }
	  wz = (tracks->dz2(itrack) <= precision2) ? precision2 : tracks->dz2(itrack);

	  double distx = pow(tx - old_x, 2) / (wx + err_x);
	  double disty = pow(ty - old_y, 2) / (wy + err_y);
	  double distz = pow(tz - old_z, 2) / (wz + err_z);

#ifdef DEBUG
	  printf("gpu fitter niter,x,dx2,wx,z,dz2,wz %d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f \n",niter,tx,wx,wx,tz,tracks->dz2(itrack),wz);
#endif

	  xpull = 0;
	  if (distz < mu2 && distx < mu2 && disty < mu2) {
	    xpull = 1.;
	    vertices->ndof(ivertex) += 1; 
	    track_weights[track_weight_counter] = xpull;
	  } else {
	    track_weights[track_weight_counter] = 0;
	  }
	  track_weight_counter++;
	  ////

	  wx = xpull / wx;
	  wy = xpull / wy;
	  wz = xpull / wz;
	  
	  x += tx * wx;
	  y += ty * wy;
	  z += tz * wz;
	  
	  s_wx += wx;
	  s_wy += wy;
	  s_wz += wz;
	  
	  s2_wx += wx * xpull;
	  s2_wy += wy * xpull;
	  s2_wz += wz * xpull;
	}
      }

      if (algorithm.useBeamConstraint) {

	wx = beamspot.cxx <=  precision2 ? 1./precision2 : 1. / (beamspot.cxx);
	wy = beamspot.cyy <=  precision2 ? 1./precision2 : 1. / (beamspot.cyy);

	x += beamspot.x * wx;
	y += beamspot.y * wy;

	s_wx  += wx;
	s2_wx += wx;
	s_wy  += wy;
	s2_wy += wy;

      }

      x /= s_wx;
      y /= s_wy;
      z /= s_wz;

      err_x = (s2_wx / (s_wx*s_wx));
      err_y = (s2_wy / (s_wy*s_wy));
      err_z = (s2_wz / (s_wz*s_wz));

      if ( (abs(x - old_x) < precision) && (abs(y - old_y) < precision) && (abs(z - old_z) < precision)){
	break;
      }
    }

    vertices->x(ivertex) = x; 
    vertices->y(ivertex) = y;
    vertices->z(ivertex) = z;

    if (algorithm.useBeamConstraint) {
      vertices->errx(ivertex) = err_x * pow(corr_x_bs,2);
      vertices->erry(ivertex) = err_y * pow(corr_x_bs,2);
      vertices->errz(ivertex) = err_z * pow(corr_z,2);
    } else {
      vertices->errx(ivertex) = err_x * pow(corr_x,2);
      vertices->erry(ivertex) = err_y * pow(corr_x,2);
      vertices->errz(ivertex) = err_z * pow(corr_z,2);
    }
    vertices->track_weight(ivertex) = track_weights;

    //----------------------------------chi2 loop-------------------------------------
    double dist = 0;
    for (unsigned int kk = 0; kk < tracks->nTrueTracks; kk++){
      unsigned int itrack = tracks->order(kk);
      if (not(tracks->isGood(itrack))) continue;

      unsigned int ivtxFromTk = tracks->kmin(itrack);
      if (ivtxFromTk == ivertex) {
	
	double wz = (tracks->dz2(itrack) <= precision) ? precision : tracks->dz2(itrack);

	if (algorithm.useBeamConstraint) {

	  double wx = (tracks->dxy2(itrack) <= precision) ? precision : tracks->dxy2(itrack);
	  double wy = (tracks->dxy2(itrack) <= precision) ? precision : tracks->dxy2(itrack);

	  dist =        pow(tracks->x(itrack) - vertices->x(ivertex), 2) /
	                (wx + vertices->errx(ivertex));
	  dist +=       pow(tracks->y(itrack) - vertices->y(ivertex), 2) /
	                (wy + vertices->erry(ivertex));

	} else {

	  double wx = (tracks->dxy2AtIP(itrack) <= precision) ? precision : tracks->dxy2AtIP(itrack);
	  double wy = (tracks->dxy2AtIP(itrack) <= precision) ? precision : tracks->dxy2AtIP(itrack);

	  dist =        pow(tracks->xAtIP(itrack) - vertices->x(ivertex), 2) /
	                (wx + vertices->errx(ivertex));
	  dist +=       pow(tracks->yAtIP(itrack) - vertices->y(ivertex), 2) /
	                (wy + vertices->erry(ivertex));
	}

	dist +=       pow(tracks->z(itrack) - vertices->z(ivertex), 2) /
                      (wz + vertices->errz(ivertex));

        chi2 += dist;
      }
    }

    vertices->chi2(ivertex) = chi2; 

#ifdef DEBUG
    printf("end of gpu fitter x,y,z %.10f,%.10f,%.10f \n",vertices->x(ivertex),vertices->y(ivertex),vertices->z(ivertex));
#endif

  }
}


#ifdef __CUDACC__
void wrapper(
    unsigned int ntracks,
    TrackForPV::TrackForPVSoA* GPUtracksObject,
    TrackForPV::VertexForPVSoA* GPUverticesObject,
    algo algorithm,
    bs beamspot
){

    //defines grid
    unsigned int blockSize = 256; //optimal size depends, probably 1 block, multiple threads per vertex
    unsigned int gridSize  = 1; //might need experimental determination

#ifdef DEBUG
    blockSize = 1;
#endif

    //action!
    fitterKernel<<<gridSize, blockSize>>>(
    	ntracks,
        GPUtracksObject,
    	GPUverticesObject,
    	algorithm,
	beamspot
    );

    cudaCheck(cudaGetLastError());
}
#endif
}




