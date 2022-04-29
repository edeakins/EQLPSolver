#include <vector>

class OCPartition{
    private:
    public:
        int target;
        int level;
        int nsplits;
        int ncsplits;
        int nrsplits;
        std::vector<int> front;
        std::vector<int> label;
        std::vector<int> unlabel;
        std::vector<int> parent;
        std::vector<int> len;
};