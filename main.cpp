#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <cstdint>
#include "lognum.h"
#include "ConversionEngine.cpp"
#define TESTCASES 10000
#define LARGEFRAC 0.1
#define BASE 2
// if you change BASE here still need to change "log2" ConversionEngine.cpp (lognum.h doesn't care, I think...)

using namespace std;


string sign(bool b);
int convertToInt(lognum L);
double convertToDouble(lognum L);
void runMSETestforDeltaComparison();
void printDemoCases();
void MACdemo();
void runTests();
void runTestsBigNums();
void runTestsSmallNums();
void runTestsMixedNums();
void runMSETest();
void runMSE_Delta_Test_Procedure();
void printMSE_Errors(int largeErrorCt,double MSE_MAC, double MSE_ACC, double MSE_MUL);
void printMSE_Errors_Short(int largeErrorCt, double MSE);

template<typename T>
void printErrors(T addErrors, T addSum, T multErrors, T multSum, int largeCt);
void printMinMax();

int main() {
    srand(time(NULL));

    runMSETest();
    runTestsSmallNums();
    runMSE_Delta_Test_Procedure();

/*
    printDemoCases();
    runTestsMixedNums();
    runTests();
    runTestsSmallNums();
    // deprecated? //runTestsBigNums();

    MACdemo();
    runMSETest();
    // deprecated DON'T USE runMSETestforDeltaComparison()
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
    double d = pow(BASE,L.getLogval());
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
    double d = pow(BASE,L.getLogval());
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

    // These tests squeeze x into range -1 < x_real < +1
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

/*
 * Functional but deprecated
 * */
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

    cout << "MSE calculated on 10k MAC operations with floating point inputs on (-2^15, 2^15)\n";

    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {
        // get integer bases, range up to 2^15
        int base_num = 1024*32; // changed for base = sqrt(2)
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

    printMSE_Errors(largeErrorCt,MSE,MSE_ACC_REF,MSE_MUL);
}

/*
 * This function was deprecated on 9/21/20
 * Do not use it.  For explanation see email to pabeerel@usc.edu on 9/23
 * which was a correction to the email of 9/22
 * */
void runMSETestforDeltaComparison() {
    cout << "STOP.\n";
    cout << "Entering a deprecated function.  Stop and find corrected one" << endl;
    int unexpectedErrorCt(0);
    double P_MSE_MUL(0.0),P_MSE(0.0),P_MSE_ACC_REF(0.0); int P_largeErrorCt(0);
    double M_MSE_MUL(0.0),M_MSE(0.0),M_MSE_ACC_REF(0.0); int M_largeErrorCt(0);

    cout << "Begin delta function testing\n";
    printf("MSE calculated on %d MAC operations with floating point inputs on (-2^15, 2^15)\n\n", (TESTCASES));

    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {

/*      // not used so we can emphasize cases where d is small and errors are worse

        // get integer bases, range up to 2^15
        int base_num = 1024*32;
        int32_t r1_int((rand() % (base_num)));
        int32_t r2_int((rand() % (base_num)));
        int32_t r3_int((rand() % (base_num)));
*/

        // get floating point values for A + B*C
        double r1 = (double)rand()/RAND_MAX; double r2 = (double)rand()/RAND_MAX; double r3 = (double)rand()/RAND_MAX;
        double pr1(r1),pr2(r2),pr3(r3),mr1(r1),mr2(r2),mr3(r3);

        // make half of the signbit pairs (pos,pos) and other half (neg,neg) for delta plus
        if (i%2) {pr2*= -1.0; pr3*= -1.0;}

        // make all of the signbit pairs (neg, pos) for delta minus
        mr2 *= -1.0;

        // log forms
        lognum PR1(toLogNum(pr1)), PR2(toLogNum(pr2)),PR3(toLogNum(pr3));
        double d_deltaplus = abs(PR2.getLogval()-PR3.getLogval());
        lognum MR1(toLogNum(mr1)), MR2(toLogNum(mr2)),MR3(toLogNum(mr3));
        double d_deltaminus = abs(MR2.getLogval()-MR3.getLogval());

        // take a peak, make sure conversions gone right
        double refPR1(convertToDouble(PR1)), refPR2(convertToDouble(PR2)),refPR3(convertToDouble(PR3)),
            refMR1(convertToDouble(MR1)), refMR2(convertToDouble(MR2)), refMR3(convertToDouble(MR3));

        //Dplus tests
        lognum p_mac_result = PR1;
        p_mac_result.MAC(PR2,PR3);

        // MSE on multiplication, addition and overall MAC
        double p_calcMUL = convertToDouble(lognum::multiplyReals(PR2,PR3));
        double p_refMUL = pr2*pr3;
        double p_MULdiff = p_refMUL-p_calcMUL;
        P_MSE_MUL += p_MULdiff*p_MULdiff;

        // this addition is not performed inside of the since this is B+C, but it is a useful comparison
        double p_acc_ref_diff = (convertToDouble(lognum::addReals(PR2,PR3))-(pr2+pr3));
        P_MSE_ACC_REF = p_acc_ref_diff * p_acc_ref_diff;

        double p_calcMAC = convertToDouble(p_mac_result);
        double p_expectedMAC = pr1 + pr2*pr3;

        double p_diff_MAC = (p_calcMAC-p_expectedMAC);
        P_MSE += p_diff_MAC*p_diff_MAC;

        if (abs(p_diff_MAC) > LARGEFRAC * abs(p_expectedMAC)) {
            ++P_largeErrorCt;
        }

        // end Dplus tests

        // Dminus tests
        lognum m_mac_result = MR1;
        m_mac_result.MAC(MR2,MR3);

        // see where errors came from
        double m_calcMUL = convertToDouble(lognum::multiplyReals(MR2,MR3));
        double m_refMUL = mr2*mr3;

        double m_MULdiff = m_refMUL-m_calcMUL;
        M_MSE_MUL += m_MULdiff*m_MULdiff;

        double m_acc_ref_diff = (convertToDouble(lognum::addReals(MR2,MR3))-(mr2+mr3));
        M_MSE_ACC_REF = m_acc_ref_diff * m_acc_ref_diff;

        double m_calcMAC = convertToDouble(m_mac_result);
        double m_expectedMAC = mr1 + mr2*mr3;

        double m_diff_MAC = (m_calcMAC-m_expectedMAC);
        M_MSE += m_diff_MAC*m_diff_MAC;

        if (abs(m_diff_MAC) > LARGEFRAC * abs(m_expectedMAC)) {
            ++M_largeErrorCt;
        }

        // end Dminus tests

        if (abs(p_acc_ref_diff) > abs(m_acc_ref_diff)) {
//            printf("Error when the inputs were: b=%f, c=%f\n",(pr2),(pr3));
            ++unexpectedErrorCt;
        }
        printf("\nDelta_Plus d value is: %f",d_deltaplus);
        printf("\nDelta_Minus d value is : %f\n",d_deltaminus);
    }

/*
    printf("Testing specifically on delta plus function - B+C operands must have same sign\n");
    printMSE_Errors(P_largeErrorCt,P_MSE,P_MSE_ACC_REF,P_MSE_MUL);

    printf("\n\nTesting specifically on delta minus function - B+C operands must have opposite sign\n");
    printMSE_Errors(M_largeErrorCt,M_MSE,M_MSE_ACC_REF,M_MSE_MUL);

    printf("\n\nUnexpected errors (d+ error > d- error) occured on %2.4f percent of cases",(100.0*unexpectedErrorCt/TESTCASES));
*/

    cout << "\nExiting a deprecated function.  Stop and find corrected one" << endl;
}

/*
 * This functions was created to replace
 * "runMSETestforDeltaComparison"
 *
 * Corresponds to "Corrections_to_Testing_Problem",
 * which was a fix to the incorrect "Testing_Problem"
 * */

void runMSE_Delta_Test_Procedure() {
    double DP_MSE(0.0),DM_MSE(0.0); int DP_bigErrorCt(0),DM_bigErrorCt(0),badErrorCt(0);
    double DP_totalPer_Error(0.0),DM_totalPer_Error(0.0);

    for (int i = 0; i < TESTCASES; ++i) {
        int base_num = 1024*32;
        int32_t r1_int((rand() % (base_num)));



        // 1. Pick a number |a|......OKAY to change a by adding an integer, making it negative, etc. (+=r1_int)
        double a = (double)rand()/RAND_MAX; //a += r1_int; leaving this part commented squeezes a_real into (-1,+1)
        // 2. Calculate A_fixed (float here)
        double Afloat = log(a)/log(BASE); // log2(a);
        // 3. Calculate B_fixed (float here) by moving up to 1 away
        double eps = (double)rand()/RAND_MAX;
        double Bfloat = Afloat + eps;
        // 4. Pick a number |b|
        double b = pow(BASE,Bfloat);

        // Exercise all test cases keeping a same for both, splitting b
        double bplus(b),bminus(b);
        if (i % 2) {
            a *= -1.0;
        }

        if (a<0) {
           if (bplus>0) {bplus *= -1.0;}
           if (bminus<0) {bminus *= -1.0;}
        } else { // a > 0 so make b > 0, bminus < 0
            if (bplus<0) {bplus *= -1.0;}
            if (bminus>0) {bminus *= -1.0;}
        }

        // what we should get
        double trueDPlus_addn = (a + bplus);
        double trueDMinus_addn = (a + bminus);

        // move to log domain & compute
        lognum A(toLogNum(a)),Bplus(toLogNum(bplus)),Bminus(toLogNum(bminus));
        lognum plusResult = lognum::addReals(A,Bplus);
        lognum minusResult = lognum::addReals(A,Bminus);

        // what we got
        double calcDPlus_addn = convertToDouble(plusResult);
        double calcDMinus_addn = convertToDouble(minusResult);

        // residual
        double DPlus_diff = abs(trueDPlus_addn - calcDPlus_addn);
        double DMinus_diff = abs(trueDMinus_addn - calcDMinus_addn);

        // percent error
        double DPlus_addn_Error_Per = 100.0 * abs(trueDPlus_addn-calcDPlus_addn)/abs(trueDPlus_addn);
        double DMinus_addn_Error_Per = 100.0 * abs(trueDMinus_addn-calcDMinus_addn)/abs(trueDMinus_addn);
        if (DPlus_addn_Error_Per > DMinus_addn_Error_Per) {++badErrorCt;}
        DP_totalPer_Error += DPlus_addn_Error_Per;
        DM_totalPer_Error += DMinus_addn_Error_Per;


        // MSE
        DP_MSE += (DPlus_diff) * (DPlus_diff);
        DM_MSE += (DMinus_diff) * (DMinus_diff);
    }
    printf("In %d cases, this # had larger delta_minus addition error: %2.6f",TESTCASES,(100.0*(TESTCASES-badErrorCt)/TESTCASES));
    cout <<"%\n";

    cout << "\nPrinting results of delta plus tests from V2 procedure" << endl;
//    printMSE_Errors_Short(DP_bigErrorCt,DP_MSE);
    printf("Avg percent error %2.6f",(DP_totalPer_Error/TESTCASES));

    cout << "%\n\nPrinting results of delta minus tests from V2 procedure" << endl;
//    printMSE_Errors_Short(DM_bigErrorCt,DM_MSE);
    printf("Avg percent error %2.6f",(DM_totalPer_Error/TESTCASES));cout <<"%";

}

void printMSE_Errors(int largeErrorCt,double MSE_MAC, double MSE_ACC, double MSE_MUL) {
    double RMSD_MAC = sqrt(MSE_MAC/TESTCASES);
    double RMSD_ACC = sqrt(MSE_ACC/TESTCASES);
    double RMSD_MUL = sqrt(MSE_MUL/TESTCASES);
    printf("Test cases off by > 10%: \t\t\t %2.6f",((100.0*largeErrorCt)/(TESTCASES))); cout << " %\n";
    printf("MSE on %d MAC operatons: \t\t\t%.4f",(TESTCASES),(MSE_MAC));
    printf("\nMSE from reference ADD ops was : \t\t%.4f",(MSE_ACC));
    printf("\nMSE from the MUL ops was : \t\t\t%.4f",(MSE_MUL));

    printf("\n\nRMSD on %d MAC operations: \t\t\t%.7f",(TESTCASES), (RMSD_MAC));
    printf("\nRMSD of ADD ops operations: \t\t\t%.7f",(RMSD_ACC));
    printf("\nRMSD of MUL operations: \t\t\t%.7f",(RMSD_MUL));cout<<endl;
}

void printMSE_Errors_Short(int largeErrorCt, double MSE) {
    double RMSD = sqrt(MSE/TESTCASES);
    printf("Test cases off by > 10%: \t %2.9f",((100.0*largeErrorCt)/(TESTCASES))); cout << " %\n";
    printf("MSE on %d ops: \t\t%.4f",(TESTCASES),(MSE));
    printf("\nRMSD on %d ops was: \t\t%.7f",(TESTCASES),(RMSD));
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