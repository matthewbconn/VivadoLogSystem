//
// Created by matth on 9/7/2020.
//

#ifndef VIVADOLOGSYSTEM_LOGNUM_H
#define VIVADOLOGSYSTEM_LOGNUM_H

#define WBITS 12
#define IBITS 6
#define NBITS 1

#include "ap_fixed.h"
//#include "hls_math.h"

typedef ap_fixed<WBITS,IBITS,AP_RND_ZERO,AP_SAT,NBITS> fixedtype;

class lognum {
public:
    // Default ctor, copy ctor, component-wise ctor
    lognum();
    lognum(const lognum &original);
    lognum(bool sign, fixedtype);

    // basic operations
    static lognum addReals(lognum v1, lognum v2);
    static lognum multiplyReals (lognum v1, lognum v2);
    void MAC(lognum A, lognum B);

    // access to components
    bool getSignBit() const;
    double getLogval();

private:
    // defining fields
    bool signBit;
    fixedtype logval;

};

// Delta function operating on a fixed type
fixedtype deltaPlus(fixedtype);
fixedtype deltaMinus(fixedtype);

#endif //VIVADOLOGSYSTEM_LOGNUM_H
