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

using Vector512d = Eigen::Matrix<double, 1024, 1>;

namespace fitterCUDA {

__global__ void fitterKernel(
    unsigned int ntracks,
    TrackForPV::TrackForPVSoA* tracks,
    TrackForPV::VertexForPVSoA* vertices,
    algo algorithm,
    bs beamspot
){

  /*
	OUTPUTS:
		DONE nTrueVertex (filled in clusterizer)
		DONE isGood flag (filled in clusterizer, can be modified (for now no mod))
		DONE x, y, z, t
		DONE chi2 (dep. on errs)
		DONE ndof (+1 per track if track has any influence)
		DONE errx, erry, errz
		DONE ntracks (filled in clusterizer)
		DONE track_id (type Vector512d, save position here. ignoring...)
		DONE track_weight (type Vector512d, save weight in order here. ignoring...)
  */


  size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x;
  size_t gridSize = blockDim.x * gridDim.x;
  double precision = 1e-24; //FLOAT CHANGE
  double precision2 = precision*precision;

  //1 block for each vertex (for block in BLOCKS)
  for (unsigned int k = firstElement; k < vertices->nTrueVertex(0); k += gridSize) {
    unsigned int ivertex = vertices->order(k);
    if (!vertices->isGood(ivertex)) continue; //skip if not good

    /*
		WEIGHTEDMEANFITTER ORDER

	fill err(x,y,z) with 4,4,400, corr(x,y,z) with 1.2,1.2,1.4
	x = sum(xi * (1/dxi^2)) / sum((1/dxi^2)), continue with y,z (check for precision)
	define err_x = 1 / sum((1/dxi^2)), continue with y,z
	for every point calc weight (xi - x) / (dxi^2 + err_x) < 9 ?
		set (ntracks, ndof) if > 0
	x = sum(xi * weight / (dxi^2)) / sum(dxi^2)
	err(x,y,z) = sum(dxi^2 * weight) * corr(x,y,z)^2 / (dxi^2)^2
	chi2 = sum((xi-x)^2/(dxi^2+err(x)) + y,z)

    */


    //ASSUMING THROUGHOUT THAT tracks->x == p.first.x() and tracks->dx == p.second.x()


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
      if not(tracks->isGood(itrack)) continue;
      unsigned int ivtxFromTk = tracks->kmin(itrack);
      if (ivtxFromTk == k) {
        //this is a valid track that is associated with this vertex
        track_ids[track_id_counter] = itrack;
        track_id_counter++;

        wx = (tracks->dxy2(itrack) <= precision2) ? 1./(precision2) : 1./tracks->dxy2(itrack);
        wy = (tracks->dxy2(itrack) <= precision2) ? 1./(precision2) : 1./tracks->dxy2(itrack);
        wz = (1./tracks->dz2(itrack) <= precision2) ? 1./(precision2) : tracks->dz2(itrack);

        x += tracks->x(itrack) * wx; 
        y += tracks->y(itrack) * wy;
        z += tracks->z(itrack) * wz;

        s_wx += wx;
        s_wy += wy;
        s_wz += wz;

	printf("gpu fitter x,dx2,wx,z,dz2,wz %f,%f,%f,%f,%f,%f \n",tracks->x(itrack),tracks->dxy2(itrack),wx,tracks->z(itrack),(1./tracks->dz2(itrack)),wz);
 
      }
    }

    vertices->track_id(ivertex) = track_ids;

    wx = beamspot.cxx <=  precision ? 1./precision2 : 1. / (beamspot.cxx*beamspot.cxx);
    wy = beamspot.cyy <=  precision ? 1./precision2 : 1. / (beamspot.cyy*beamspot.cyy);

    x += beamspot.x * wx;
    y += beamspot.y * wy;

    x /= (s_wx + wx);
    y /= (s_wy + wy);
    z /= s_wz;

    double old_x, old_y, old_z;

    double xpull;
    int niter = 0;
    double mu = 3.;
    double mu2 = mu*mu;

    double err_x, err_y, err_z;

    err_x = 1. / s_wx;
    err_y = 1. / s_wy;
    err_z = 1. / s_wz;

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
        if not(tracks->isGood(itrack)) continue;
	unsigned int ivtxFromTk = tracks->kmin(itrack);
	if (ivtxFromTk == k) {

	  ////WEIGHTING
	  double wx = (tracks->dxy2(itrack) <= precision2) ? precision2 : tracks->dxy2(itrack);
	  double wy = (tracks->dxy2(itrack) <= precision2) ? precision2 : tracks->dxy2(itrack);
	  double wz = (1./tracks->dz2(itrack) <= precision2) ? precision2 : 1./tracks->dz2(itrack);
	  
	  double distx = pow(tracks->x(itrack) - old_x, 2) / (wx + err_x);
	  double disty = pow(tracks->y(itrack) - old_y, 2) / (wy + err_y);
	  double distz = pow(tracks->z(itrack) - old_z, 2) / (wz + err_z);

	  xpull = 0;
	  if (distz < mu2 && distx < mu2 && disty < mu2) {
	    xpull = 1.;
	    vertices->ndof(ivertex) += 1; 
	    //vertices->ntracks(ivertex) += 1; already given by clusterizer
	    track_weights[track_weight_counter] = xpull;
	    printf("fitterCUDA: adding track with z=%f to vtx. %d \n",tracks->z(itrack),k);
	  } else {
	    track_weights[track_weight_counter] = 0;
	  }
	  track_weight_counter++;
	  ////

	  wx = xpull / wx;
	  wy = xpull / wy;

	  wz = xpull / wz;
	  
	  x += tracks->x(itrack) * wx;
	  y += tracks->y(itrack) * wy;
	  z += tracks->z(itrack) * wz;
	  
	  s_wx += wx;
	  s_wy += wy;
	  s_wz += wz;
	  
	  s2_wx += wx * xpull;
	  s2_wy += wy * xpull;
	  s2_wz += wz * xpull;
	}
      }

      // WHY IS IT DIFFERENT THAN THE OTHER TIME AROUND???
      wx = beamspot.cxx <=  precision2 ? 1./precision2 : 1. / (beamspot.cxx);
      wy = beamspot.cyy <=  precision2 ? 1./precision2 : 1. / (beamspot.cyy);

      x += beamspot.x * wx;
      y += beamspot.y * wy;

      s_wx  += wx;
      s2_wx += wx;
      s_wy  += wy;
      s2_wy += wy;

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

    vertices->errx(ivertex) = err_x * pow(corr_x_bs,2);
    vertices->erry(ivertex) = err_y * pow(corr_x_bs,2);
    vertices->errz(ivertex) = err_z * pow(corr_z,2);

    vertices->track_weight(ivertex) = track_weights;

    //----------------------------------chi2 loop-------------------------------------
    //creates error: GammaContinuedFraction::a too large, ITMAX too small
    //maybe just average chi2s from the tracks? those are being calculated anyway.

    double dist = 0;
    for (unsigned int kk = 0; kk < tracks->nTrueTracks; kk++){
      unsigned int itrack = tracks->order(kk);
      if not(tracks->isGood(itrack)) continue;
      unsigned int ivtxFromTk = tracks->kmin(itrack);
      if (ivtxFromTk == k) {

        double wx = (tracks->dxy2(itrack) <= precision) ? precision : tracks->dxy2(itrack);
        double wy = (tracks->dxy2(itrack) <= precision) ? precision : tracks->dxy2(itrack);
        double wz = (1./tracks->dz2(itrack) <= precision) ? precision : 1./tracks->dz2(itrack);

        //printf("wx,wy,wz: %f, %f, %f\n", wx, wy, wz);
        //printf("errx,erry,errz: %f, %f, %f\n", vertices->errx(ivertex), vertices->erry(ivertex), vertices->errz(ivertex));

        dist =        pow(tracks->x(itrack) - vertices->x(ivertex), 2) /
                      (wx + vertices->errx(ivertex));
        dist +=       pow(tracks->y(itrack) - vertices->y(ivertex), 2) /
                      (wy + vertices->erry(ivertex));
        dist +=       pow(tracks->z(itrack) - vertices->z(ivertex), 2) /
                      (wz + vertices->errz(ivertex));
        //printf("%f\n", dist);
        chi2 += dist;
      }
    }

    vertices->chi2(ivertex) = chi2; //requirement
    //weight tracks by chi2?

    printf("end of fitter: %d, %f, %f, %f, %d, %f, %f\n",
      ivertex,
      vertices->x(ivertex),
      vertices->y(ivertex),
      vertices->z(ivertex),
      vertices->ntracks(ivertex),
      vertices->errz(ivertex),
      vertices->chi2(ivertex)
    );

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
    unsigned int blockSize = 1; //optimal size depends, probably 1 block, multiple threads per vertex
    unsigned int gridSize  = 1; //might need experimental determination

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




