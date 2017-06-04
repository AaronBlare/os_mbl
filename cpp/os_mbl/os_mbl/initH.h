#ifndef __INIT_H__
#define __INIT_H__

#include "Model.h"

void init_hE_vector(Model * m, ConfigData &cd, ConfigParam &cp);
void init_hU_vector(Model * m, ConfigData &cd, ConfigParam &cp);
void init_hJ_vector(Model * m, ConfigData &cd, ConfigParam &cp);

void init_H0(Model * m);

#endif