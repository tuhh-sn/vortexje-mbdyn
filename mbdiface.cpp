#include "mbdiface.h"
using namespace std;

mbdInterface::mbdInterface() {

    SocketPath = "127.0.0.1";
        port = 7771;
        printf("Socket: %s\n",SocketPath);

        verbose = 0;
        flg_dbg = 0;

        /** Connect(2) timeout.  When peer is not listening:
         * - a value of 0 results in immediate return with error;
         * - a positive value corresponds to a timeout in seconds;
         * - a negative value corresponds to never timing out.
         */
        timeout = -1;
        data_and_next = 1;

        rigid = 1;
        labels = 0;
        rot = 512;
        accels = 0;
}

mbc_nodal_t* mbdInterface::Init() {
        /* assign mbc handles */
        mbcx = { { 0 } };           // empty nodal struct
        _mbcHandle = &mbcx;         // pointer to nodal struct

        _mbcHandle->mbc.timeout = timeout;
        _mbcHandle->mbc.verbose = verbose;
        _mbcHandle->mbc.data_and_next = data_and_next;

        cout << "Connecting socket..." << endl;
        socketflag = mbc_inet_init(reinterpret_cast<mbc_t*>(_mbcHandle), SocketPath, port); // const_cast for function argument passing

        if(socketflag == 0){
            printf("Socket successfully connected!\n");
        };

        socketflag = mbc_nodal_init(_mbcHandle, rigid, nodes, labels, rot, accels);
        cout << "Nodal init flag = " << socketflag << endl;

        socketflag = mbc_nodal_negotiate_request(_mbcHandle);
        cout << "Nodal negotiation = " << socketflag << endl;

        return _mbcHandle;
}


mbc_nodal_t* mbdInterface::getMotion(void) {
    socketflag = mbc_nodal_get_motion(_mbcHandle);

    if (socketflag<0) {
    #ifdef MEX_COMPILE_FLAG
        mexErrMsgTxt("getMotion failed.");
        exit;
    #endif
    }

    return _mbcHandle;
}


int mbdInterface::setForce(double *_real, unsigned last) {
    if (flg_dbg) {
        cout << "Settings forces..." << endl;
    }

    for (int k=0; k<nodes; k++) {
        for (int j=0; j<3; j++) {
            if (flg_dbg) {
                cout << "F: input[" << j << "," << k << "]: " << _real[j+k*6] << endl;
                cout << "M: input[" << j << "," << k << "]: " << _real[j+k*6+3] << endl;
            }

            /* write forces */
            MBC_N_F(_mbcHandle)[j+k*3] = _real[j+k*6];
            if (flg_dbg) {
                cout << "F: output[" << j+k*3 << "]: " << _mbcHandle->n_d_f[j+k*3] << endl;
            }

            /* write moments */
            MBC_N_M(_mbcHandle)[j+k*3] = _real[j+k*6+3];
            if (flg_dbg) {
                cout << "M: output[" << j+k*3 << "]: " << _mbcHandle->n_d_m[j+k*3] << endl;
            }
        }
    }


    mbc_nodal_put_forces(_mbcHandle, last);
    if (flg_dbg) {
        cout << "Nodal setForce = " << socketflag << endl;
    }
    if (socketflag<0) {
        #ifdef MEX_COMPILE_FLAG
            mexErrMsgTxt("setForce failed.");
        #endif
        exit;
    }

    return socketflag;
}

/* Eigen vector input */
int mbdInterface::setForce(Eigen::Vector3d &forces, Eigen::Vector3d &moments, unsigned last) {
    Eigen::VectorXd real(6);
    real.head(3) = forces;
    real.tail(3) = moments;
    return mbdInterface::setForce(real.data(), last);
}

mbdInterface::~mbdInterface() {
    cout << "Closing mbdInterface." << endl;
    mbc_nodal_destroy(_mbcHandle);
}
