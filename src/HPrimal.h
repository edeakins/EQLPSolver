#ifndef HPRIMAL_H_
#define HPRIMAL_H_

#include "HModel.h"

class HPrimal {
public:
    void solvePhase2(HModel *ptr_model);
private:
    void primalRebuild();
    void primalChooseColumn();
    void primalChooseRow();
    void primalUpdate();
    void primalChooseLinker();
    void primalCollect();

    // Model pointer
    HModel *model;
    int numCol;
    int numRow;
    int numTot;

    // Pivot related
    int limitUpdate;
    int countUpdate;
    int invertHint;
    int columnIn;
    int rowOut;

    // Solve buffer
    HVector row_ep;
    HVector row_ap;
    HVector column;
    double row_epDensity;
    double columnDensity;
};

#endif /* HPRIMAL_H_ */
