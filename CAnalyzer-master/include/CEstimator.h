#ifndef C_ESTIMATOR_H
#define C_ESTIMATOR_H

#include <vector>
#include "CMatrix.h"

class TF1;
class TFormula;

class CEstimator
{
public:
    struct DataPoint
    {
        double x;
        double val;
        double error;

        DataPoint()
        : x(0.), val(0.), error(0.)
        {};

        DataPoint(const double &xi, const double &vi, const double &err)
        : x(xi), val(vi), error(err)
        {};
    };

    struct Parameter
    {
        double initial;
        double prev_value;
        double value;
        double step;
        double base_step;
        bool lock;

        Parameter()
        : initial(0.), prev_value(0.), value(0.),
          step(0.), base_step(0.), lock(true)
        {};
        Parameter(const double &v, const double &s = 0.)
        : initial(v), prev_value(v), value(v), step(0), lock(true)
        {
            if(s == 0.)
                base_step = value*0.01;
        }
    };

public:
    CEstimator(const std::string &path = "");
    virtual ~CEstimator();

    // manipulate members
    void LoadFormula(const std::string &path);
    void SaveFormula(const std::string &path);
    void DeleteFormula();
    void UnlockPar(const size_t &i) {parameters[i].lock = false;};
    void LockPar(const size_t &i) {parameters[i].lock = true;};
    void UpdatePars();
    void SetFormula(TF1 *tf);
    void SetFormula(const char *c);
    void SetDataPoints(const std::vector<double> &x,
                       const std::vector<double> &y,
                       const std::vector<double> &err);
    void SetDataPoints(const size_t &n,
                       const double *x,
                       const double *y,
                       const double *err);
    void GenerateCovMatrix();
    void SetParameter(const size_t &i,
                      const double &p,
                      const double &step = 0.);
    void SetParameters(const std::vector<double> &p);
    void SetStepFactor(const double &fine, const double &coarse);

    // Maximum Likelihood:
    // weight matrix <- covariance matrix
    // penalty matrix <- comes from physical meaning to restrict the change

    // Weighted Least Square
    // weight matrix <- weighting the data points
    // penalty matrix <- none

    // Least Square
    // weight matrix <- unity matrix
    // penalty matrix <- none
    void SetWeightMatrix(const CMatrix &m);
    void SetPenaltyMatrix(const CMatrix &m);

    // get members
    size_t GetNpar() const {return parameters.size();};
    size_t GetNpoints() const {return data.size();};
    TFormula *GetFormula();
    const std::vector<DataPoint> &GetDataPoints() const {return data;};
    const std::vector<Parameter> &GetParameters() const {return parameters;};
    Parameter GetParameter(const size_t &i) const
    {
        if(i >= parameters.size())
            return Parameter(0);
        return parameters.at(i);
    }
    double GetFormulaVal(const double &x);

    // fit related
    // fit, iteration limitation, step range
    virtual void Fit(int c_iter = 100, int f_iter = 100, int range = 30, bool verbose = true);
    virtual bool Optimize(int range, double factor, bool verbose);
    virtual double Evaluate(const double &factor = 0.);
    virtual void NextStep(const double &factor, bool verbose);
    virtual CMatrix GetHessian();
    virtual void CalcStep();

    // fit quality check
    virtual double GetReducedChiSquare();
    virtual double GetPearsonChiSquare();
    virtual double GetAbsoluteError();
    virtual double GetNLL_Gaussian();
    virtual double GetRootMeanSquaredError();
    virtual double GetAkaikeCriterion();
    virtual double GetBayesianCriterion();
    virtual double GetHannanCriterion(const double &c);
    virtual double GetKashyapCriterion(const CMatrix &F_M);

private:
    TFormula *formula;
    CMatrix M_weight;
    bool weight_diagonal;
    CMatrix M_penalty;
    bool penalty_diagonal;
    double fine_step_size;
    double coarse_step_size;
    std::vector<DataPoint> data;
    std::vector<Parameter> parameters;
};

#endif
