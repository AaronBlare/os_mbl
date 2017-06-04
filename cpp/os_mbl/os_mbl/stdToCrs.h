#ifndef __STD_TO_CRS__
#define __STD_TO_CRS__

#include "Model.h"
#include <vector>

using namespace std;

crsMatrix * stdToCrs(vector<pair<int , dcomplex> > * mat, int N);

#endif