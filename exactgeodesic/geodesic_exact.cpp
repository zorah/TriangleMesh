/**
 *
 * Compile with: mex -v geodesic_exact.cpp
 * Indices are expected to be in cpp format (starting at 0).
 *
 */

#include "mex.h"
#include "geodesic_algorithm_exact.h"
#include "geodesic_error.h"
#include "geodesic_datastructures.h"
#include <iostream>

#include <cstring>
#include <string>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
  // DEBUG
	ofstream fout("./debug.txt", ios::out);
	fout << nlhs << " " << nrhs << "\n";

// check and read mex input
  if (nrhs != 4 || nlhs != 2)
    mexErrMsgTxt("Usage: [distances, sources] = geodesic_exact(vertices, triangles, sources, targets). Transform indices in triangles to initial 0.");

  // read vertices
  const double* const vertices = mxGetPr(prhs[0]);
  const long nVertices = long( mxGetN(prhs[0]) );
  const int dimVertices = int( mxGetM(prhs[0]) );

  // read triangles
  const double* const triangles = mxGetPr(prhs[1]);
  const long nTriangles = long( mxGetN(prhs[1]) );
  const int dimTriangles = int( mxGetM(prhs[1]) );

  // read sources
  const double* const sources = mxGetPr(prhs[2]);
  const int nSources = int( mxGetN(prhs[2]) );
  const int dimSources = int( mxGetM(prhs[2]) );
  
  // read targets
  const double* const targets = mxGetPr(prhs[3]);
  const int nTargets = int( mxGetN(prhs[3]) );
  const int dimTargets = int( mxGetM(prhs[3]) );

  // check dimensions
  if (dimTriangles != 3)
    mexErrMsgTxt("Triangles must be given in 3xm.");

  if (dimVertices != 3)
    mexErrMsgTxt("Vertices must be given in 3xn.");
  
  if (dimSources != 1)
    mexErrMsgTxt("Sources must be given in 1xk.");
  
  if (dimTargets != 1)
    mexErrMsgTxt("Targets must be given in 1xk.");

    fout << "Check complete." << "\n";

// calculate geodesics
try
{
  geodesic::Mesh mesh_;
  geodesic::GeodesicAlgorithmExact* algorithm_;
  geodesicPoints* points_ = new geodesicPoints(vertices, int(nVertices));
  geodesicFaces* faces_ = new geodesicFaces(triangles, int(nTriangles));

  fout << "Declaration complete." << "\n";

  mesh_.initialize_mesh_data(*points_, *faces_);

  fout << "Mesh initialization complete." << "\n";

  algorithm_ = new geodesic::GeodesicAlgorithmExact(&mesh_);

  fout << "Algorithm initialization complete." << "\n";

  geodesic::SurfacePoint platzhalter(&mesh_.vertices()[ 0 ]);
  std::vector<geodesic::SurfacePoint> all_sources(nSources, platzhalter);
  for (int i=0; i < nSources; i++) {
      all_sources[i] = geodesic::SurfacePoint(&mesh_.vertices()[ sources[i] ]);
  }
  algorithm_->propagate(all_sources);

  fout << "Propagation complete." << "\n";

  // read out result
  plhs[0] = mxCreateDoubleMatrix(nTargets, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(nTargets, 1, mxREAL);
  double* distmap = mxGetPr(plhs[0]);
  double* bestSources = mxGetPr(plhs[1]);

  for(int i=0; i < nTargets; i++)
  {
    geodesic::SurfacePoint p(&mesh_.vertices()[targets[i] ]);
    double distance;
    bestSources[i] = double(algorithm_->best_source(p, distance));
    distmap[i] = distance;
  }
} catch (std::exception& e)
  {
    const std::string msg = "Exception caught: " + std::string(e.what());
    mexErrMsgTxt(msg.c_str());
  }

}
