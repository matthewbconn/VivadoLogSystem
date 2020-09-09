//
// Created by matth on 9/7/2020.
//

#include "lognum.h"

// with no better information, generate for x(real) = 1
// so sgn(x) = 1, log2(x) = 0
lognum::lognum(): signBit(true), logval(0) {}

// ctor to instantiate for a specific value (someone else's job to convert)
lognum::lognum(bool sign, fixedtype fixedNum):
        signBit(sign), logval(fixedNum){}

lognum lognum::addReals(lognum v1, lognum v2) {
    lognum result;

    result.signBit = v1.logval > v2.logval ? v1.signBit : v2.signBit;

    fixedtype d = v1.logval > v2.logval ?
            v1.logval - v2.logval : v2.logval - v1.logval;

    fixedtype larger = v1.logval > v2.logval ? v1.logval : v2.logval;

    result.logval  = v1.signBit == v2.signBit ?
            (larger + deltaPlus(d)): (larger + deltaMinus(d));

    return result;
}

lognum lognum::multiplyReals(lognum v1, lognum v2) {
    return lognum(! v1.signBit xor v2.signBit, v1.logval + v2.logval);
}

/*
 * Red flag: you coded this when you were tired.
 *           that's probably why it isn't working/
 * */
lognum lognum::MAC(lognum A, lognum B) {
    lognum multiply(multiplyReals(A,B));
    lognum result(addReals(*this,multiply));
    return (addReals(*this,multiplyReals(A,B)));
}

bool lognum::getSignBit() const  {
    return signBit;
}

double lognum::getLogval() {
    return logval.to_double();
    // this is returning an integer, not a floating point value
    return double(logval);
}

fixedtype deltaPlus(fixedtype d) {
    fixedtype bitshift(1);
    bitshift = (bitshift >> d);
    return bitshift;
}

fixedtype deltaMinus(fixedtype d) {
    fixedtype bitshift(1.5);
    bitshift = -1*(bitshift >> d);
    return bitshift;
}
