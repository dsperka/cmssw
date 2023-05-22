#ifndef GeometryVector_GlobalPoint_h
#define GeometryVector_GlobalPoint_h

#include "DataFormats/GeometryVector/interface/GlobalTag.h"
#include "DataFormats/GeometryVector/interface/Point3DBase.h"

typedef Point3DBase<float, GlobalTag> Global3DPoint;
typedef Point3DBase<double, GlobalTag> Global3DPointDouble;

// Global points are three-dimensional by default
typedef Global3DPoint GlobalPoint;
typedef Global3DPointDouble GlobalPointDouble;

#endif  // GeometryVector_GlobalPoint_h
