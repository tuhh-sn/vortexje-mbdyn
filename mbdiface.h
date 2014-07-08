#ifndef MBDINTERFACE_H
#define MBDINTERFACE_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <mbc.h>
#include <assert.h>
#include <sstream>
#include <math.h>
#include <eigen3/Eigen/Eigen>

class mbdInterface {
private:
    const char *SocketPath;
    short unsigned port, verbose, timeout, data_and_next, rigid, labels, rot, accels;
    bool flg_dbg;
    int socketflag;
    mbc_nodal_t	mbcx;                           // nodal handler
    mbc_nodal_t * _mbcHandle;                   // nodal pointer

public:
    mbdInterface();             // constructor
    mbc_nodal_t* Init();        // Init interface
    short unsigned nodes;
    mbc_nodal_t* getMotion(void);
    int setForce(double *_real, unsigned last);
    void setForce(double real, unsigned last);
    int setForce(Eigen::Vector3d &forces, Eigen::Vector3d &moments, unsigned last);
    ~mbdInterface();            // destructor
};

#endif // MBDINTERFACE_H
