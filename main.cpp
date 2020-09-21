#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <cstdint>
#include "lognum.h"
#include "ConversionEngine.cpp"
#define TESTCASES 10000
#define LARGEFRAC 0.1

using namespace std;


string sign(bool b);
int convertToInt(lognum L);
double convertToDouble(lognum L);
void printDemoCases();
void MACdemo();
void runTests();
void runTestsBigNums();
void runTestsSmallNums();
void runTestsMixedNums();
void runMSETest();

template<typename T>
void printErrors(T addErrors, T addSum, T multErrors, T multSum, int largeCt);
void printMinMax();

int main() {
    srand(time(NULL));
    runMSETest();
/*
    printDemoCases();
    runTestsMixedNums();
    runTests();
    runTestsSmallNums();
    // deprecated? //runTestsBigNums();

    MACdemo();

    printMinMax();
*/

    return 0;
}

string sign(bool b) {
    if(b) {return ("positive");}
    return("negative");
}

/*
 * Input: a lognum object
 * Output: the corresponding real value, as an int
 * */
int convertToInt(lognum L) {
    double temp = L.getLogval();
    double d = pow(2,L.getLogval());
    if (L.getSignBit()) {
        return d;
    }
    return (int)(-1*d);
}

/*
 * Input: a lognum object
 * Output: the corresponding real value, as a double
 * */
double convertToDouble(lognum L) {
    double d = pow(2,L.getLogval());
    if (L.getSignBit()) {
        return d;
    }
    return (-1.0*d);
}

void printDemoCases() {
    int v(0), w(1), x(4),y(-2),z(-8);
    cout << "x = " << x << endl;
    cout << "y = " << y << endl;

    lognum V(toLogNum(v)), W(toLogNum(w)),X(toLogNum(x)),Y(toLogNum(y)),Z(toLogNum(z));
    lognum XplusY = lognum::addReals(X,Y);
    lognum XtimesY = lognum::multiplyReals(X,Y);

    cout << "X + Y is: " << sign(XplusY.getSignBit()) << endl;
    cout << "X + Y logval: " << XplusY.getLogval() << endl;
    cout << "X + Y realval: " << convertToDouble(XplusY) << endl;
    cout << "X + Y expected: " << x+y << "\n" << endl;


    cout << "X * Y is: " << sign(XtimesY.getSignBit()) << endl;
    cout << "X * Y logval: " << XtimesY.getLogval() << endl;
    cout << "X * Y realval: " << convertToDouble(XtimesY) << endl;
    cout << "X * Y expected: " << x*y << "\n" << endl;

    lognum ZplusY(lognum::addReals(Z,W));
    lognum ZtimesY(lognum::multiplyReals(Z,W));
    cout << "v = " << v << endl;
    cout << "w = " << w << endl;
    cout << "z = " << z << endl;
    cout << "Z + W to Integer computed as: " << convertToInt(ZplusY) << endl;
    cout << "Z + W to Float computed as: " << convertToDouble(ZplusY) << endl;
    cout << "Z * W to Integer computed as: " << convertToInt(ZtimesY) << endl;
    cout << "Z * W to Float computed as: " << convertToDouble(ZtimesY) << endl << endl;

    cout << "0 + 1 as integer: " << convertToInt(lognum::addReals(V,W)) << endl;
    cout << "0 + 1 as double: " << convertToDouble(lognum::addReals(V,W)) << endl;
    cout << "0 * 1 as integer: " << convertToInt(lognum::multiplyReals(V,W)) << endl;
    cout << "0 * 1 as double: " << convertToDouble(lognum::multiplyReals(V,W)) << endl;
}

void MACdemo() {
    int m(1); lognum M(toLogNum(m));
    int n(5); lognum N(toLogNum(n));
    cout << "Initial value: " << convertToDouble(N) << "\n";
    int numits(5);
    for (int i = 0; i < numits; ++i) {
        N.MAC(M,M);
    }
    cout << "After " << numits << " multiply accumulate ops: " << convertToDouble(N) << "\n" << endl;
}

void runTests() {
    double addErrors(0), addErrorMargin(0), multErrors(0), multErrorMargin(0);
    int largeMulErrorCt(0);

    cout << "runTests\n";

    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {
        int32_t r1((rand() % (1024*32)));
        int32_t r2((rand() % (1024*32)));
        lognum R1(toLogNum(r1)), R2(toLogNum(r2));

        // ensure conversion back and forth works right
        int r1test(convertToInt(R1)),r2test(convertToInt(R2));

        long long expSum = r2 + r1;
        double calcSum = convertToDouble(lognum::addReals(R1,R2));
        addErrors += abs(calcSum - 1.0*expSum);
        double addMargin = abs(expSum) < abs(calcSum) ? abs(expSum) : abs(calcSum);
        addErrorMargin += addMargin;

        double expProd = 1.0*r2 * r1;
        double calcProd(convertToDouble(lognum::multiplyReals(R1,R2)));

        multErrors += abs(calcProd - 1.0*expProd);
        if (abs(calcProd - expProd) > LARGEFRAC * abs(expProd)) {
            ++largeMulErrorCt;
        }
        double multMargin = abs(expProd) < abs(calcProd) ? abs(expProd) : abs(calcProd);
        multErrorMargin += multMargin;
    }

    cout << "Resuts for both numbers big: " << endl;
    printErrors(addErrors,addErrorMargin,multErrors,multErrorMargin,largeMulErrorCt);
}

void runTestsSmallNums() {
    double addErrors(0), addErrorMargin(0), multErrors(0), multErrorMargin(0);
    int largeMulErrorCt(0);

    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {
        double r1(1.0*((rand() % 1000))/1000); if (i%2) {r1 *= -1;}
        double r2(1.0*((rand() % 1000))/1000); if (i%4) {r2 *= -1;}
        lognum R1(toLogNum(r1)), R2(toLogNum(r2));

        // ensure conversion back and forth works fine
        double r1test(convertToDouble(R1)),r2test(convertToDouble(R2));
        double expSum = r2 + r1;
        double expProd = r2 * r1;
        double calcSum = convertToDouble(lognum::addReals(R1,R2));
        double calcProd = convertToDouble(lognum::multiplyReals(R1,R2));
        addErrors += abs(calcSum - 1.0*expSum);
        double addMargin = abs(expSum) < abs(calcSum) ? abs(expSum) : abs(calcSum);
        addErrorMargin += addMargin;
        multErrors += abs(calcProd - 1.0*expProd);
        if (abs(calcProd - 1.0 * expProd) > LARGEFRAC * abs(expProd)) {
            ++largeMulErrorCt;
        }
        double multMargin = abs(expProd) < abs(calcProd) ? abs(expProd) : abs(calcProd);
        multErrorMargin += multMargin;
    }

    cout << "Results for both numbers small:" << endl;
    printErrors(addErrors,addErrorMargin,multErrors,multErrorMargin,largeMulErrorCt);
}

void runTestsMixedNums() {
    double addErrors(0), addErrorMargin(0), multErrors(0), multErrorMargin(0);
    int largeMulErrorCt(0);

    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {
        double r1(1.0*((rand() % 1000))/1000); if (i%2) {r1 *= -1;}
        // change the parens
        double r2((rand() % 1024*1024) + 0);
        lognum R1(toLogNum(r1)), R2(toLogNum(r2));

        // ensure conversion back and forth works fine
        double r1test(convertToDouble(R1)),r2test(convertToDouble(R2));
        double expSum = r2 + r1;
        double expProd = r2 * r1;
        double calcSum = convertToDouble(lognum::addReals(R1,R2));
        double calcProd = convertToDouble(lognum::multiplyReals(R1,R2));
        addErrors += abs(calcSum - 1.0*expSum);
        double addMargin = abs(expSum) < abs(calcSum) ? abs(expSum) : abs(calcSum);
        addErrorMargin += addMargin;
        multErrors += abs(calcProd - 1.0*expProd);
        if (abs(calcProd - 1.0 * expProd) > LARGEFRAC * abs(expProd)) {
            ++largeMulErrorCt;
        }
        double multMargin = abs(expProd) < abs(calcProd) ? abs(expProd) : abs(calcProd);
        multErrorMargin += multMargin;
    }

    cout << "Results for one number small one number big:" << endl;
    printErrors(addErrors,addErrorMargin,multErrors,multErrorMargin,largeMulErrorCt);
}

void runTestsBigNums() {
    double addErrors(0), addErrorMargin(0), multErrors(0), multErrorMargin(0);
    int largeMulErrorCt(0);

    cout << "runTestsBigNums\n";
    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {
        int r1((rand() % 1024*32));
        int r2((rand() % 1024*32));
        lognum R1(toLogNum(r1)), R2(toLogNum(r2));

        // ensure conversion back and forth works right
        int r1test(convertToInt(R1)),r2test(convertToInt(R2));

        int expSum = r2 + r1;
        double calcSum = convertToDouble(lognum::addReals(R1,R2));
        addErrors += abs(calcSum - 1.0*expSum);
        double addMargin = abs(expSum) < abs(calcSum) ? abs(expSum) : abs(calcSum);
        addErrorMargin += addMargin;

        double expProd = 1.0*r2 * r1;

        // why is calcProd not changing at all?! calcSum changes
        double calcProd = convertToDouble(lognum::multiplyReals(R1,R2));

        multErrors += abs(calcProd - 1.0*expProd);
        if (abs(calcProd - expProd) > LARGEFRAC * abs(expProd)) {
            ++largeMulErrorCt;
        }
        double multMargin = abs(expProd) < abs(calcProd) ? abs(expProd) : abs(calcProd);
        multErrorMargin += multMargin;
    }

    cout << "Resuts for both numbers big: " << endl;
    printErrors(addErrors,addErrorMargin,multErrors,multErrorMargin,largeMulErrorCt);
}

void runMSETest() {
    double MSE_MUL(0.0),MSE(0.0),MSE_ACC_REF(0.0); int largeErrorCt(0);

    cout << "MSE calculated on 10k MAC operations with floating point inputs on (-2^15, 2^15)\n\n";

    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {
        // get integer bases, range up to 2^15
        int base_num = 1024*32;
        int32_t r1_int((rand() % (base_num)));
        int32_t r2_int((rand() % (base_num)));
        int32_t r3_int((rand() % (base_num)));

        // get floating point
        double r1 = static_cast<double>(r1_int) + (double)rand()/RAND_MAX;
        double r2 = static_cast<double>(r2_int) + (double)rand()/RAND_MAX;
        double r3 = static_cast<double>(r3_int) + (double)rand()/RAND_MAX;

        // make some of them negative
        if (i%2) {r1*= -1.0;}
        if (i%3) {r2*= -1.0;}
        if (i%4) {r3*= -1.0;}

        lognum R1(toLogNum(r1)), R2(toLogNum(r2)),R3(toLogNum(r3));

        // ensure conversion back and forth works right
        int r1test(convertToDouble(R1)),r2test(convertToDouble(R2)),r3test(convertToDouble(R3));

        lognum mac_res = R1;
        mac_res.MAC(R2,R3);

        // see where errors came from
        double calcBase = convertToDouble(R1);
        double calcMUL = convertToDouble(lognum::multiplyReals(R2,R3));
        double refMUL = r2*r3;
        double MULdiff = refMUL-calcMUL;
        double acc_ref_diff = (convertToDouble(lognum::addReals(R2,R3))-(r2+r3));

        MSE_MUL += MULdiff*MULdiff;
        MSE_ACC_REF = acc_ref_diff * acc_ref_diff;

        double calcMAC = convertToDouble(mac_res);
        double expectedMAC = r1 + r2*r3;

        double diff = (calcMAC-expectedMAC);
        MSE += diff*diff;

        if (abs(diff) > LARGEFRAC * abs(expectedMAC)) {
            ++largeErrorCt;
        }
    }

    double RMSD = sqrt(MSE/TESTCASES);
    double RMSD_MUL = sqrt(MSE_MUL/TESTCASES);
    double RMSD_ACC_REF = sqrt(MSE_ACC_REF/TESTCASES);

    printf("Test cases off by > 10%: \t\t\t %2.6f",((100.0*largeErrorCt)/(TESTCASES))); cout << " %\n";
    printf("MSE on 10k MAC operatons: \t\t\t%.2f",(MSE));
    printf("\nMSE from reference ADD ops was : \t\t%.2f",(MSE_ACC_REF));
    printf("\nMSE from the MUL ops was : \t\t\t%.2f",(MSE_MUL));

    printf("\n\nRMSD on 10k MAC operations: \t\t\t%.5f",(RMSD));
    printf("\nRMSD of 10k reference ADD ops operations: \t%.5f",(RMSD_ACC_REF));
    printf("\nRMSD of 10k MUL operations: \t\t\t%.5f",(RMSD_MUL));
}

template<typename T>
void printErrors(T addErrors, T addSum, T multErrors, T multSum, int largeCt) {
    cout << "Avg. error per addition test\t" << addErrors/TESTCASES << endl;
    printf("Relative to range\t\t %2.9f",((100.0*addErrors)/(addSum*TESTCASES))); cout << " %\n";
    cout << "Avg. error per mult. test\t" << multErrors/TESTCASES << endl;
    printf("Relative to range\t\t %2.9f",((100.0*multErrors)/(multSum*TESTCASES))); cout << " %\n";
    printf("Percent large errors\t\t %2.9f",((100.0*largeCt)/(TESTCASES))); cout << " %\n\n";
}

void printMinMax() {
    cout << "Min is: " << convertToDouble(toLogNum(INT64_MIN)) << endl;
    cout << "Max is: " << convertToDouble(toLogNum(INT64_MAX)) << endl;
}