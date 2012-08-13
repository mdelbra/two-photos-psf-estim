//Copyright (C) 2007 Lionel Moisan: initial version
//Copyright (C) 2010-2011 Pascal Monasse, Pierre Moulon: genericity, C++ class
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "libOrsa/orsa_model.hpp"
#include "libOrsa/conditioning.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace orsa {

/// Points are normalized according to image dimensions. Subclass must
/// initialize logalpha0_ and take account of normalization in doing that.
OrsaModel::OrsaModel(const Mat &x1, int w1, int h1,
                     const Mat &x2, int w2, int h2)
: x1_(x1.nrow(), x1.ncol()), x2_(x2.nrow(), x2.ncol()),
  N1_(3,3), N2_(3,3), logalpha0_(0) {
  assert(2 == x1_.nrow());
  assert(x1_.nrow() == x2_.nrow());
  assert(x1_.ncol() == x2_.ncol());

  NormalizePoints(x1, &x1_, &N1_, w1, h1);
  NormalizePoints(x2, &x2_, &N2_, w2, h2);
}

/// If multiple solutions are possible, return false.
bool OrsaModel::ComputeModel(const std::vector<size_t> &indices,
                             Model *model) const {
    std::vector<Model> models;
    Fit(indices, &models);
    if(models.size() != 1)
      return false;
    *model = models.front();
    Unnormalize(model);
    return true;
}

/// logarithm (base 10) of binomial coefficient
static float logcombi(int k, int n)
{
  if (k>=n || k<=0) return(0.0);
  if (n-k<k) k=n-k;
  double r = 0.0;
  for (int i = 1; i <= k; i++)
    r += log10((double)(n-i+1))-log10((double)i);

  return static_cast<float>(r);
}

/// tabulate logcombi(.,n)
static void makelogcombi_n(int n, std::vector<float> & l)
{
  l.resize(n+1);
  for (int k = 0; k <= n; k++)
    l[k] = logcombi(k,n);
}

/// tabulate logcombi(k,.)
static void makelogcombi_k(int k,int nmax, std::vector<float> & l)
{
  l.resize(nmax+1);
  for (int n = 0; n <= nmax; n++)
    l[n] = logcombi(k,n);
}

/// Find best NFA and number of inliers wrt square error threshold in e.
OrsaModel::ErrorIndex OrsaModel::bestNFA(const std::vector<ErrorIndex>& e,
                                         double loge0,
                                         double maxThreshold,
                                         const std::vector<float> &logc_n,
                                         const std::vector<float> &logc_k) const
{
  const int startIndex = SizeSample();
  const double multError = (DistToPoint()? 1.0: 0.5);

  ErrorIndex bestIndex(std::numeric_limits<double>::infinity(), startIndex);
  const size_t n = e.size();
  for(size_t k=startIndex+1; k<=n && e[k-1].first<=maxThreshold; ++k) {
    float logalpha = logalpha0_ + multError *log10(e[k-1].first);
    ErrorIndex index(loge0 +
                     logalpha * (double)(k-startIndex) +
                     logc_n[k] +
                     logc_k[k], k);

    if(index.first < bestIndex.first)
      bestIndex = index;
  }
  return bestIndex;
}

/// Get a (sorted) random sample of size X in [0:n-1]
static void random_sample(std::vector<size_t> &k, int X, size_t n)
{
  for(size_t i=0; i < (size_t)X; i++) {
    size_t r = (rand()>>3)%(n-i), j;
    for(j=0; j<i && r>=k[j]; j++)
      r++;
    size_t j0 = j;
    for(j=i; j > j0; j--)
      k[j]=k[j-1];
    k[j0] = r;
  }
}

/// Pick a random sample
/// \param sizeSample The size of the sample.
/// \param vec_index  The possible data indices.
/// \param sample The random sample of sizeSample indices (output).
static void UniformSample(int sizeSample,
                          const std::vector<size_t> &vec_index,
                          std::vector<size_t> *sample) {
  sample->resize(sizeSample);
  random_sample(*sample, sizeSample, vec_index.size());
  for(int i = 0; i < sizeSample; ++i)
    (*sample)[i] = vec_index[ (*sample)[i] ];
}

/// Generic implementation of 'ORSA':
/// A Probabilistic Criterion to Detect Rigid Point Matches
///    Between Two Images and Estimate the Fundamental Matrix.
/// Bibtex :
/// @article{DBLP:journals/ijcv/MoisanS04,
///  author    = {Lionel Moisan and B{\'e}renger Stival},
///  title     = {A Probabilistic Criterion to Detect Rigid Point Matches
///    Between Two Images and Estimate the Fundamental Matrix},
///  journal   = {International Journal of Computer Vision},
///  volume    = {57},
///  number    = {3},
///  year      = {2004},
///  pages     = {201-218},
///  ee        = {http://dx.doi.org/10.1023/B:VISI.0000013094.38752.54},
///  bibsource = {DBLP, http://dblp.uni-trier.de}
///}
/// 
/// ORSA is based on an a contrario criterion of
/// inlier/outlier discrimination, is parameter free and relies on an optimized
/// random sampling procedure. It returns the log of NFA and optionally
/// the best estimated model.
///
/// \param vec_inliers Output vector of inlier indices.
/// \param nIter The number of iterations.
/// \param precision (input/output) threshold for inlier discrimintation.
/// \param model The best computed model.
/// \param bVerbose Display optimization statistics.
double OrsaModel::orsa(std::vector<size_t> & vec_inliers,
                       size_t nIter,
                       double *precision,
                       Model *model,
                       bool bVerbose) const {
  vec_inliers.clear();

  const int sizeSample = SizeSample();
  const size_t nData = x1_.ncol();
  if(nData <= (size_t)sizeSample)
    return std::numeric_limits<double>::infinity();

  const double maxThreshold = (precision && *precision>0)?
    *precision * *precision *N2_(0,0)*N2_(0,0): // Square max error
    std::numeric_limits<double>::infinity();

  std::vector<ErrorIndex> vec_residuals(nData); // [residual,index]
  std::vector<size_t> vec_sample(sizeSample); // Sample indices

  // Possible sampling indices (could change in the optimization phase)
  std::vector<size_t> vec_index(nData);
  for (size_t i = 0; i < nData; ++i)
    vec_index[i] = i;

  // Precompute log combi
  double loge0 = log10((double)NbModels() * (nData-sizeSample));
  std::vector<float> vec_logc_n, vec_logc_k;
  makelogcombi_n(nData, vec_logc_n);
  makelogcombi_k(sizeSample,nData, vec_logc_k);

  // Output parameters
  double minNFA = std::numeric_limits<double>::infinity();
  double errorMax = 0;

  // Main estimation loop.
  for (size_t iter=0; iter < nIter; iter++) {
    UniformSample(sizeSample, vec_index, &vec_sample); // Get random sample

    std::vector<Model> vec_models; // Up to max_models solutions
    Fit(vec_sample, &vec_models);

    // Evaluate models
    for (size_t k = 0; k < vec_models.size(); ++k)
    {
      // Residuals computation and ordering
      for (size_t i = 0; i < nData; ++i)
      {
        double error = Error(vec_models[k], i);
        vec_residuals[i] = ErrorIndex(error, i);
      }
      std::sort(vec_residuals.begin(), vec_residuals.end());

      // Most meaningful discrimination inliers/outliers
      ErrorIndex best = bestNFA(vec_residuals, loge0, maxThreshold,
                                vec_logc_n, vec_logc_k);

      if(best.first < 0 && best.first < minNFA) // A better model was found
      {
        minNFA = best.first;
        vec_inliers.resize(best.second);
        for (size_t i=0; i<best.second; ++i)
          vec_inliers[i] = vec_residuals[i].second;
        errorMax = vec_residuals[best.second-1].first; // Error threshold
        if(model) *model = vec_models[k];

        // ORSA optimization: draw samples among best set of inliers so far
        vec_index = vec_inliers;
        if(bVerbose)
        {
          std::cout << "  nfa=" << minNFA
                    << " inliers=" << vec_index.size()
                    << " precision=" << sqrt(errorMax)/N2_(0,0)
                    << " (iter=" << iter;
          std::cout << ",sample=" << vec_sample.front();
          std::vector<size_t>::const_iterator it=vec_sample.begin();
          for(++it; it != vec_sample.end(); ++it)
            std::cout << ',' << *it;
          std::cout << ")" <<std::endl;
        }
      }
    }
  }
  if(precision)
    *precision = sqrt(errorMax)/N2_(0,0);
  if(model && !vec_inliers.empty())
    Unnormalize(model);
  return minNFA;
}

} // namespace orsa
