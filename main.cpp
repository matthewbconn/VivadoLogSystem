#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include "lognum.h"
#include "ConversionEngine.cpp"
#define TESTCASES 10000

using namespace std;


string sign(bool b);
int convertToInt(lognum L);
double convertToDouble(lognum L);
void printDemoCases();
void MACdemo();
void runTests();
void runTestsSmallNums();
void runTestsMixedNums();

int main() {
    cout << "Hello, World!" << "\n" << endl;

    printDemoCases();
    MACdemo();

    srand(time(NULL));
    runTests();
    runTestsSmallNums();
    runTestsMixedNums();

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
    int d = pow(2,L.getLogval());
    if (L.getSignBit()) {
        return d;
    }
    return (-1*d);
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
    int n(0); lognum N(toLogNum(n));
    for (int i = 0; i < 10; ++i) {
        N.MAC(M,M);
    }
    cout << "\nAfter 10 multiply accumulate ops: " << convertToDouble(N) << "\n" << endl;
}

void runTests() {
    double addErrors(0);
    double addErrorMargin(0);
    double multErrors(0);
    double multErrorMargin(0);
    int largeMulErrorCt(0);

    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {
        int r1((rand() % 1024*1024) - 127);
        int r2((rand() % 1024*1024) + 0);
        lognum R1(toLogNum(r1)), R2(toLogNum(r2));
        // notice this line right here: you are getting the next lowest power of 2...wonder why
        int r1test(convertToInt(R1)),r2test(convertToInt(R2));
        int expSum = r2 + r1;
        int expProd = r2 * r1;
        double calcSum = convertToDouble(lognum::addReals(R1,R2));
        double calcProd = convertToDouble(lognum::multiplyReals(R1,R2));
        addErrors += abs(calcSum - 1.0*expSum);
        double addMargin = abs(expSum) < abs(calcSum) ? abs(expSum) : abs(calcSum);
        addErrorMargin += addMargin;
        multErrors += abs(calcProd - 1.0*expProd);
        if (abs(calcProd - 1.0 * expProd) > 0.1 * abs(expProd)) {
//            cout << "Large error when multiplying " << r1 << " to " << r2 << endl;
            ++largeMulErrorCt;
        }
        double multMargin = abs(expProd) < abs(calcProd) ? abs(expProd) : abs(calcProd);
        multErrorMargin += multMargin;
    }

    cout << "Resuts for both numbers big: " << endl;
    cout << "Avg. error per addition test\t" << addErrors/TESTCASES << endl;
    printf("Relative to range\t\t %2.4f",((100.0*addErrors)/(addErrorMargin*TESTCASES))); cout << " %\n";
    cout << "Avg. error per mult. test\t" << multErrors/TESTCASES << endl;
    printf("Relative to range\t\t %2.4f",((100.0*multErrors)/(multErrorMargin*TESTCASES))); cout << " %\n";
    cout << "Total large errors: " << largeMulErrorCt << "\n\n";
}

void runTestsSmallNums() {
    double addErrors(0);
    double addErrorMargin(0);
    double multErrors(0);
    double multErrorMargin(0);
    int largeMulErrorCt(0);

    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {
        double r1(1.0*((rand() % 1000))/1000); if (i%2) {r1 *= -1;}
        double r2(1.0*((rand() % 1000))/1000); if (i%4) {r2 *= -1;}
        lognum R1(toLogNum(r1)), R2(toLogNum(r2));
        // notice this line right here: you are getting the next lowest power of 2...wonder why
        double r1test(convertToDouble(R1)),r2test(convertToDouble(R2));
        double expSum = r2 + r1;
        double expProd = r2 * r1;
        double calcSum = convertToDouble(lognum::addReals(R1,R2));
        double calcProd = convertToDouble(lognum::multiplyReals(R1,R2));
        addErrors += abs(calcSum - 1.0*expSum);
        double addMargin = abs(expSum) < abs(calcSum) ? abs(expSum) : abs(calcSum);
        addErrorMargin += addMargin;
        multErrors += abs(calcProd - 1.0*expProd);
        if (abs(calcProd - 1.0 * expProd) > 0.1 * abs(expProd)) {
//            cout << "Large error when multiplying " << r1 << " to " << r2 << endl;
            ++largeMulErrorCt;
        }
        double multMargin = abs(expProd) < abs(calcProd) ? abs(expProd) : abs(calcProd);
        multErrorMargin += multMargin;
    }

    cout << "Results for both numbers small:" << endl;
    cout << "Avg. error per addition test\t" << addErrors/TESTCASES << endl;
    printf("Relative to range\t\t %2.9f",((100.0*addErrors)/(addErrorMargin*TESTCASES))); cout << " %\n";
    cout << "Avg. error per mult test\t" << multErrors/TESTCASES << endl;
    printf("Relative to range\t\t %2.9f",((100.0*multErrors)/(multErrorMargin*TESTCASES))); cout << " %\n";
    cout << "Total large errors: " << largeMulErrorCt << "\n\n";
}

void runTestsMixedNums() {
    double addErrors(0);
    double addErrorMargin(0);
    double multErrors(0);
    double multErrorMargin(0);
    int largeMulErrorCt(0);

    // Currently allotting 6 signed bits for integer
    // log2(|x|) can range from 2^-5 ( = -32)  to just below 2^5 - 1 ( = 31 - eps)
    // Corresponding x ranges from  just above 0 to just below 2^31 ~ 10^9
    for (int i = 0; i < TESTCASES; ++i) {
        double r1(1.0*((rand() % 1000))/1000); if (i%2) {r1 *= -1;}
        double r2((rand() % 1024*1024) + 0);
        lognum R1(toLogNum(r1)), R2(toLogNum(r2));
        // notice this line right here: you are getting the next lowest power of 2...wonder why
        double r1test(convertToDouble(R1)),r2test(convertToDouble(R2));
        double expSum = r2 + r1;
        double expProd = r2 * r1;
        double calcSum = convertToDouble(lognum::addReals(R1,R2));
        double calcProd = convertToDouble(lognum::multiplyReals(R1,R2));
        addErrors += abs(calcSum - 1.0*expSum);
        double addMargin = abs(expSum) < abs(calcSum) ? abs(expSum) : abs(calcSum);
        addErrorMargin += addMargin;
        multErrors += abs(calcProd - 1.0*expProd);
        if (abs(calcProd - 1.0 * expProd) > 0.1 * abs(expProd)) {
//            cout << "Large error when multiplying " << r1 << " to " << r2 << endl;
            ++largeMulErrorCt;
        }
        double multMargin = abs(expProd) < abs(calcProd) ? abs(expProd) : abs(calcProd);
        multErrorMargin += multMargin;
    }

    cout << "Results for one number small one number big:" << endl;
    cout << "Avg. error per addition test\t" << addErrors/TESTCASES << endl;
    printf("Relative to range\t\t %2.4f",((100.0*addErrors)/(addErrorMargin*TESTCASES))); cout << " %\n";
    cout << "Avg. error per mult test\t" << multErrors/TESTCASES << endl;
    printf("Relative to range\t\t %2.4f",((100.0*multErrors)/(multErrorMargin*TESTCASES))); cout << " %\n";
    cout << "Total large errors: " << largeMulErrorCt << "\n\n";
}