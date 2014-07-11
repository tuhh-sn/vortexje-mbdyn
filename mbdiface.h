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
    int socketflag;
    mbc_nodal_t	mbcx;                           // nodal handler
    mbc_nodal_t * _mbcHandle;                   // nodal pointer

public:
    mbdInterface();
    mbc_nodal_t* Init();        // Init interface

    /**
     * @brief nodes
     *
     * number of nodes, each has 6 dof's
     */
    short unsigned nodes;

    /**
     * @brief getMotion update MBDyn buffer of kinematics
     * @return mbc_nodal_t handle
     */
    mbc_nodal_t* getMotion(void);

    /**
     * @brief setForce write forces into MBDyn buffer and send convergence flag
     * @param _real pointer to double array of length nodes
     * @param last convergence flag
     * @return 0 on success
     */
    int setForce(double *_real, unsigned last);

    /**
     * @brief setForce
     * @param real double array
     * @param last
     */
    void setForce(double real, unsigned last);

    /**
     * @brief setForce
     * @param forces Eigen 3d vector of forces
     * @param moments ...and moments
     * @param last
     * @return 0 on success
     */
    int setForce(Eigen::Vector3d &forces, Eigen::Vector3d &moments, unsigned last);
    ~mbdInterface();
};

#endif // MBDINTERFACE_H
