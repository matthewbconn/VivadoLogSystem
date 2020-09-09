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
    lognum();
    lognum(bool sign, fixedtype);
    static lognum addReals(lognum v1, lognum v2);
    static lognum multiplyReals (lognum v1, lognum v2);
    lognum MAC(lognum A, lognum B);

    bool getSignBit() const;
    double getLogval();

private:
    bool signBit;
    fixedtype logval;

};

fixedtype deltaPlus(fixedtype);
fixedtype deltaMinus(fixedtype);

#endif //VIVADOLOGSYSTEM_LOGNUM_H
