#include <iostream>
#include <map>
#include <vector>
#include <random>
#include <memory>
#include <stdexcept>
#include <cmath>
#include <iomanip>

const unsigned RAND_SEED = std::random_device()();
//const unsigned RAND_SEED = 4084505731;
double random_c()
{
    static std::mt19937 rand(RAND_SEED);
    return 1.0*(1 + rand()) / (rand.max() + 2ull);
}

class DES final
{
public:
    using time_t = double;
private:
    // -------------------------------------------------------------------------------------
    class Event
    {
    protected:
        static void brokeA(DES *des) { des->brokeA(); }
        static void repairA(DES *des) { des->repairA(); }

        static void brokeB(DES *des) { des->brokeB(); }
        static void repairB(DES *des) { des->repairB(); }

        static void beginRepairing(DES *des) { des->beginRepairing(); }
        static void endRepairing(DES *des) { des->endRepairing(); }

        static void addEvent(DES* des, time_t time, std::unique_ptr<Event> pEvent) { des->addEvent(time, std::move(pEvent)); }

        static time_t genTau(double lam) { return -1.0/lam*std::log(random_c()); }
    public:
        virtual void exec(DES *des) const = 0;
        virtual ~Event() {}
    };
    friend class Event;

    // -------------------------------------------------------------------------------------
    class BrokeAEvent : public Event
    {
    public:
        virtual void exec(DES* des) const override 
        { 
            if (des->isFail())
            {
                addEvent(des, des->getCurrTime() + genTau(des->getLamA()), std::make_unique<BrokeAEvent>());
                return;
            }

            brokeA(des);
            if (des->getCorrA() >= des->getNA())
                addEvent(des, des->getCurrTime() + genTau(des->getLamA()), std::make_unique<BrokeAEvent>());
            if (!des->isRepairing())
            {
                beginRepairing(des);
                if (des->getFailA() > des->getFailB())
                    addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairAEvent>());
                else if (des->getFailA() < des->getFailB())
                    addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairBEvent>());
                else
                {
                    if (des->getLamA() > des->getLamB())
                        addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairAEvent>());
                    else
                        addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairBEvent>());
                }
            }
        }
    };
    class RepairAEvent : public Event
    {
    public:
        virtual void exec(DES* des) const override 
        { 
            endRepairing(des);
            repairA(des); 
            if (des->getCorrA() <= des->getNA())
                addEvent(des, des->getCurrTime() + genTau(des->getLamA()), std::make_unique<BrokeAEvent>());
            if (des->getFailA() > 0 || des->getFailB() > 0)
            {
                beginRepairing(des);
                if (des->getFailA() > des->getFailB())
                    addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairAEvent>());
                else if (des->getFailA() < des->getFailB())
                    addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairBEvent>());
                else
                {
                    if (des->getLamA() > des->getLamB())
                        addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairAEvent>());
                    else
                        addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairBEvent>());
                }
            }
        }
    };

    class BrokeBEvent : public Event
    {
    public:
        virtual void exec(DES* des) const override 
        { 
            if (des->isFail())
            {
                addEvent(des, des->getCurrTime() + genTau(des->getLamB()), std::make_unique<BrokeBEvent>());
                return;
            }

            brokeB(des);
            if (des->getCorrB() >= des->getNB())
                addEvent(des, des->getCurrTime() + genTau(des->getLamB()), std::make_unique<BrokeBEvent>());
            if (!des->isRepairing())
            {
                beginRepairing(des);
                if (des->getFailA() > des->getFailB())
                    addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairAEvent>());
                else if (des->getFailA() < des->getFailB())
                    addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairBEvent>());
                else
                {
                    if (des->getLamA() > des->getLamB())
                        addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairAEvent>());
                    else
                        addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairBEvent>());
                }
            }
        }
    };
    class RepairBEvent : public Event
    {
    public:
        virtual void exec(DES* des) const override 
        { 
            endRepairing(des);
            repairB(des); 
            if (des->getCorrB() <= des->getNB())
                addEvent(des, des->getCurrTime() + genTau(des->getLamB()), std::make_unique<BrokeBEvent>());
            if (des->getFailA() > 0 || des->getFailB() > 0)
            {
                beginRepairing(des);
                if (des->getFailA() > des->getFailB())
                    addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairAEvent>());
                else if (des->getFailA() < des->getFailB())
                    addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairBEvent>());
                else
                {
                    if (des->getLamA() > des->getLamB())
                        addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairAEvent>());
                    else
                        addEvent(des, des->getCurrTime() + genTau(des->getLamS()), std::make_unique<RepairBEvent>());
                }
            }
        }
    };

    // -------------------------------------------------------------------------------------
    class StartEvent : public Event
    {
    public:
        virtual void exec(DES* des) const override
        {
            for (int k = 0; k < des->getNA(); ++k)
                addEvent(des, genTau(des->getLamA()), std::make_unique<BrokeAEvent>());
            for (int k = 0; k < des->getNB(); ++k)
                addEvent(des, genTau(des->getLamB()), std::make_unique<BrokeBEvent>());
        }
    };
    class FinishEvent : public Event
    {
    public:
        virtual void exec(DES* des) const override {}
    };

    // -------------------------------------------------------------------------------------
    using fec_t = std::multimap<time_t, std::unique_ptr<Event>>;
    using cec_t = std::vector<std::unique_ptr<Event>>;

    // -------------------------------------------------------------------------------------
    double lamA, lamB;
    double lamS;
    int nA, nB;
    int rA, rB;
    time_t maxTime;

    // -------------------------------------------------------------------------------------
    fec_t fec;
    time_t currTime;
    int corrA, failA;
    int corrB, failB;
    bool repairing;

    // -------------------------------------------------------------------------------------
    void brokeA()
    {
        if (corrA <= 0) throw std::runtime_error("Can't broke A, all A are fault");
        --corrA;
        ++failA;
    }
    void repairA()
    {
        if (failA <= 0) throw std::runtime_error("Can't repair A, all A are correct");
        ++corrA;
        --failA;
    }

    void brokeB()
    {
        if (corrB <= 0) throw std::runtime_error("Can't broke B, all B are fault");
        --corrB;
        ++failB;
    }
    void repairB()
    {
        if (failB <= 0) throw std::runtime_error("Can't repair B, all B are correct");
        ++corrB;
        --failB;
    }

    // -------------------------------------------------------------------------------------
    void beginRepairing() 
    {  
        if (repairing) throw std::runtime_error("Repairing is already began");
        repairing = true;
    }
    void endRepairing()
    {
        if (!repairing) throw std::runtime_error("Repairing is already not doing");
        repairing = false;
    }

    // -------------------------------------------------------------------------------------
    void addEvent(time_t time, std::unique_ptr<Event> pEvent) 
    { 
        if (time < currTime) throw std::runtime_error("Time of the event is less than current");
        fec.emplace(time, std::move(pEvent)); 
    }

    // -------------------------------------------------------------------------------------
public:
    DES(double lamA, double lamB, double lamS, int nA, int nB, int rA, int rB, time_t maxTime)
        : lamA(lamA),lamB(lamB), lamS(lamS), nA(nA), nB(nB), rA(rA), rB(rB), maxTime(maxTime),
        currTime(0.0), corrA(nA+rA), failA(0), corrB(nB+rB), failB(0), repairing(false)
    {
        if (lamA <= 0.0 || lamB <= 0.0 || lamS <= 0.0) throw std::runtime_error("Lambdas must be possitive");
        if (nA < 1 || nB < 1 || rA < 1 || rB < 1) throw std::runtime_error("Ns and Rs must not be less then 1");
        if (maxTime <= 0.0) throw std::runtime_error("Maximum time must be possitive");

        addEvent(0.0, std::make_unique<StartEvent>());
        addEvent(maxTime, std::make_unique<FinishEvent>());
    }

    // -------------------------------------------------------------------------------------
    double getLamA() const { return lamA; }
    double getLamB() const { return lamB; }
    double getLamS() const { return lamS; }

    int getNA() const { return nA; }
    int getNB() const { return nB; }
    int getRA() const { return rA; }
    int getRB() const { return rB; }

    time_t getMaxTime() const { return maxTime; }

    // -------------------------------------------------------------------------------------
    time_t getCurrTime() const { return currTime; }
    bool isFinish() const { return currTime >= maxTime; }

    int getCorrA() const { return corrA; }
    int getFailA() const { return failA; }

    int getCorrB() const { return corrB; }
    int getFailB() const { return failB; }

    bool isRepairing() const { return repairing; }

    bool isFail() const { return corrA < 1 || corrB < nB; }

    // -------------------------------------------------------------------------------------
    bool nextTime()
    {
        if (isFinish()) return false;
        if (fec.empty()) throw std::runtime_error("Finish was not reached");

        cec_t cec;
        currTime = fec.begin()->first;
        do
        {
            cec.push_back(std::move(fec.begin()->second));
            fec.erase(fec.begin());
        }
        while(!fec.empty() && fec.begin()->first == currTime);

        for (size_t k = 0; k < cec.size(); ++k)
        {
            cec[k]->exec(this);
        }
        cec.clear();

        return true;
    }
};

// -------------------------------------------------------------------------------------
struct ResultsTable
{
    std::vector<DES::time_t> time;
    std::vector<int> corrA;
    std::vector<int> corrB;
    std::vector<bool> isRepairing;
    std::vector<bool> isFail;
};

struct Statistics
{
    double aveQ;
    double aveCountA;
    double aveCountB;
    double aveRepair;
};

// -------------------------------------------------------------------------------------
void printResults(const ResultsTable &res)
{
    for (int k = 0; k < res.time.size(); ++k)
    {
        std::cout << std::setw(10) << std::setprecision(4) << std::fixed << res.time[k]
                  << std::setw(4) << res.corrA[k]
                  << std::setw(4) << res.corrB[k]
                  << ";\n";
    }
}

Statistics getStatistics(const ResultsTable& res)
{
    Statistics stats = { 0.0, 0.0, 0.0, 0.0 };
    for (int k = 0; k < res.time.size() - 1; ++k)
    {
        DES::time_t dt = res.time[k+1] - res.time[k];

        stats.aveQ += res.isFail[k] * dt;
        stats.aveCountA += res.corrA[k] * dt;
        stats.aveCountB += res.corrB[k] * dt;
        stats.aveRepair += res.isRepairing[k] * dt;
    }
    DES::time_t time = res.time.back();
    stats.aveQ /= time; 
    stats.aveCountA /= time; 
    stats.aveCountB /= time; 
    stats.aveRepair /= time;

    return stats;
}

// -------------------------------------------------------------------------------------
int main()
{
    const double LAM_A = 3.0;
    const double LAM_B = 5.0;
    const double LAM_S = 16.0;

    const int NA = 3;
    const int NB = 2;
    const int RA = 2;
    const int RB = 1;

    const DES::time_t TIME = 1.0470;

    ResultsTable res;
    try
    {
        DES des(LAM_A, LAM_B, LAM_S, NA, NB, RA, RB, TIME);
        while (des.nextTime())
        {
            res.time.push_back(des.getCurrTime());
            res.corrA.push_back(des.getCorrA());
            res.corrB.push_back(des.getCorrB());
            res.isRepairing.push_back(des.isRepairing());
            res.isFail.push_back(des.isFail());
        }
    }
    catch (std::exception &error)
    {
        std::cerr << error.what() << "\n";
        return -1;
    }
    printResults(res);

    Statistics stats = getStatistics(res);
    std::cout << "\n"
              << "aveQ = " << stats.aveQ 
              << ", aveCountA = " << stats.aveCountA 
              << ", aveCountB = " << stats.aveCountB 
              << ", aveRepair = " << stats.aveRepair 
              << "\n";

    std::cout << "\nSeed = " << RAND_SEED << "\n";

    return 0;
}