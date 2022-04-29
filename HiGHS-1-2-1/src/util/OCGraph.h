class OCGraph{
    private:
    public:
        int numCol_;
        int numRow_;
        int nnz_;
        int numWeights_;
        int numTot_;
        int* adj;
        int* edg;
        int* wght;
        int* colors;
        int* edgColors;
};
