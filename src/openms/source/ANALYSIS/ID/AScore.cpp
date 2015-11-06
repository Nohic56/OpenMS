// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar, Timo Sachsenberg, Petra Gutenbrunner $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/AScore.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <cmath>
#include <algorithm> //find
#include <boost/math/special_functions/binomial.hpp>
#include <iostream>

using namespace std;

namespace OpenMS
{

  AScore::AScore()
  {
  }

  AScore::~AScore()
  {
  }

  PeptideHit AScore::compute(PeptideHit & hit, RichPeakSpectrum & real_spectrum, double fragment_mass_tolerance, bool fragment_mass_unit_ppm) const
  {
    PeptideHit phospho = hit;
    
    //reset phospho
    phospho.setScore(0);
    
    //String with_phospho_string = phospho.getSequence().toString();
    String sequence_str = phospho.getSequence().toString();
    
    Size number_of_phosphorylation_events = numberOfPhosphoEvents_(sequence_str);
    AASequence seq_without_phospho = removePhosphositesFromSequence_(sequence_str);
  
    // determine all phospho sites
    vector<Size> sites(getSites_(seq_without_phospho));
    Size number_of_STY = sites.size();
    
    if (number_of_STY < number_of_phosphorylation_events)
    {
      number_of_phosphorylation_events = number_of_STY;
    }
    
    if (real_spectrum.empty())
    {
      return phospho;
    }

    vector<vector<Size> > permutations(computePermutations_(sites, (Int)number_of_phosphorylation_events));
    vector<RichPeakSpectrum> th_spectra;
    if (permutations.empty())
    {
      th_spectra = createTheoreticalSpectra_(seq_without_phospho);
    }
    else
    {
      th_spectra = createTheoreticalSpectra_(permutations, seq_without_phospho);
    }
        
    // prepare real spectrum windows
    if (!real_spectrum.isSorted())
    {
      real_spectrum.sortByPosition();
    }
    vector<RichPeakSpectrum> windows_top10(peakPickingPerWindowsInSpectrum_(real_spectrum));
    
    // calculate peptide score for each possible phospho site permutation
    vector<vector<double> > peptide_site_scores(calculatePermutationPeptideScores_(th_spectra, windows_top10, fragment_mass_tolerance, fragment_mass_unit_ppm));
    
    // rank peptide permutations ascending
    multimap<double, Size> ranking(rankWeightedPermutationPeptideScores_(peptide_site_scores));
    
    phospho.setScore(ranking.rbegin()->first); // initialize score with highest peptide score (aka highest weighted score)
    phospho.setSequence(AASequence::fromString(th_spectra[ranking.rbegin()->second].getName()));
    phospho.setMetaValue("Search_engine_sequence", hit.getSequence().toString());
    
    if (number_of_phosphorylation_events == 0 || number_of_STY == 0 || number_of_STY == number_of_phosphorylation_events)
    {
      return phospho;
    }
    
    vector<ProbablePhosphoSites> phospho_sites;
    determineHighestScoringPermutations_(peptide_site_scores, phospho_sites, permutations, ranking);
    
    Int rank = 1;
    for (vector<ProbablePhosphoSites>::iterator s_it = phospho_sites.begin(); s_it < phospho_sites.end(); ++s_it)
    {
      vector<RichPeakSpectrum> site_determining_ions;
      // Previously, the precursor charge was used here. This is clearly wrong and it is better to use charge 1 here.
      computeSiteDeterminingIons_(th_spectra, *s_it, site_determining_ions);
      
      Size N = site_determining_ions[0].size(); // all possibilities have the same number so take the first one
      double p = static_cast<double>(s_it->peak_depth) / 100.0;
      
      Size n_first = 0;
      for (Size depth = 0; depth != windows_top10.size(); ++depth) // for each 100 m/z window
      {
        n_first += numberOfMatchedIons_(site_determining_ions[0], windows_top10[depth], s_it->peak_depth, fragment_mass_tolerance, fragment_mass_unit_ppm);
      }

      double P_first = computeCumulativeScore_(N, n_first, p);

      Size n_second = 0;
      for (Size depth = 0; depth <  windows_top10.size(); ++depth) //each 100 m/z window
      {
        n_second += numberOfMatchedIons_(site_determining_ions[1], windows_top10[depth], s_it->peak_depth, fragment_mass_tolerance, fragment_mass_unit_ppm);
      }

      double P_second = computeCumulativeScore_(N, n_second, p);
      
      //abs is used to avoid -0 score values
      double score_first = abs(-10 * log10(P_first));
      double score_second = abs(-10 * log10(P_second));
      double AScore_first = score_first - score_second;
      
      phospho.setMetaValue("AScore_" + String(rank), AScore_first);
      ++rank;      
    }
    
    return phospho;
  }

  double AScore::computeCumulativeScore_(Size N, Size n, double p) const
  {
    OPENMS_PRECONDITION(n <= N, "The number of matched ions (n) can be at most as large as the number of trials (N).");
    OPENMS_PRECONDITION(p >= 0 && p <= 1.0, "p must be a probability [0,1].");

    // return bad p value if none has been matched (see Beausoleil et al.)
    if (n == 0) return 1.0;

    double score = 0.0;
    // score = sum_{k=n..N}(\choose{N}{k}p^k(1-p)^{N-k})
    for (Size k = n; k <= N; ++k)
    {
      double coeff = boost::math::binomial_coefficient<double>((unsigned int)N, (unsigned int)k);
      double pow1 = pow((double)p, (int)k);
      double pow2 = pow(double(1 - p), double(N - k));
      score += coeff * pow1 * pow2;
    }

    return score;
  }

  void AScore::determineHighestScoringPermutations_(const std::vector<std::vector<double> >& peptide_site_scores, std::vector<ProbablePhosphoSites>& sites, const vector<vector<Size> >& permutations, std::multimap<double, Size>& ranking) const
  {
    // For every phospho site of the highest (weighted) scoring phospho site assignment:
    // 1. determine the next best (weighted) score assignment with this site in unphosporylated state.
    // 2. determine the filtering level (peak depths) that maximizes the (unweighted) score difference between these two assignments

    sites.clear();
    // take first set of phospho site assignments
    sites.resize(permutations[0].size());    
    const vector<Size> & best_peptide_sites = permutations[ranking.rbegin()->second]; // sites of the assignment that achieved the highest weighted score

    for (Size i = 0; i < best_peptide_sites.size(); ++i)  // for each phosphorylated site
    {
      multimap<double, Size>::reverse_iterator rev = ranking.rbegin();      
      sites[i].first = best_peptide_sites[i]; // store the site
      sites[i].seq_1 = rev->second; // and permutation
      bool peptide_not_found = true;
      
      // iterate from best scoring peptide to the first peptide that doesn't contain the current phospho site
      do
      {      
        ++rev;
        for (Size j = 0; j < best_peptide_sites.size(); ++j)
        {
          if (j == i)
          {
            if (find(permutations[rev->second].begin(), permutations[rev->second].end(), best_peptide_sites[j]) != permutations[rev->second].end())
            {
              peptide_not_found = true;
              break;
            }
            else
            {
              peptide_not_found = false;
            }
          }
          else
          {
            if (find(permutations[rev->second].begin(), permutations[rev->second].end(), best_peptide_sites[j]) == permutations[rev->second].end())
            {
              peptide_not_found = true;
              break;
            }
            else
            {
              peptide_not_found = false;
            }
          }
        }
      }
      while (peptide_not_found);

      // store permutation of peptide without the phospho site i (seq_2)
      sites[i].seq_2 = rev->second;

      // store phospho site location that is not contained in the best scoring (seq_1) but in seq_2.
      for (Size j = 0; j < permutations[sites[i].seq_2].size(); ++j)
      {
        if (find(permutations[sites[i].seq_1].begin(), permutations[sites[i].seq_1].end(), permutations[sites[i].seq_2][j]) == permutations[sites[i].seq_1].end())
        {
          sites[i].second = permutations[sites[i].seq_2][j];
          break;
        }
      }
    }

    // store peak depth that achieves maximum score difference between best and runner up for every phospho site.
    for (Size i = 0; i < sites.size(); ++i)
    {
      double maximum_score_difference = 0.0;
      sites[i].peak_depth = 1;
      vector<double>::const_iterator first_it = peptide_site_scores[sites[i].seq_1].begin();
      vector<double>::const_iterator second_it = peptide_site_scores[sites[i].seq_2].begin();
      
      for (Size depth = 1; second_it < peptide_site_scores[sites[i].seq_2].end(); ++second_it, ++first_it, ++depth)
      {
        double phospho_at_site_score = *first_it;
        double no_phospho_at_site_score = *second_it;
        double score_difference = phospho_at_site_score - no_phospho_at_site_score;
        
        if (score_difference > maximum_score_difference)
        {
          maximum_score_difference = score_difference;
          sites[i].peak_depth = depth;
        }
      }
    }
  }
  
  // calculation of the number of different speaks between the theoretical spectra of the two best scoring peptide permutations, respectively
  void AScore::computeSiteDeterminingIons_(vector<RichPeakSpectrum> & th_spectra, ProbablePhosphoSites & candidates, vector<RichPeakSpectrum> & site_determining_ions) const
  {
    site_determining_ions.clear();
    site_determining_ions.resize(2);
    TheoreticalSpectrumGenerator spectrum_generator;
    
    RichPeakSpectrum spectrum_first = th_spectra[candidates.seq_1];
    RichPeakSpectrum spectrum_second = th_spectra[candidates.seq_2];
    
    RichPeakSpectrum spectrum_first_diff;
    AScore::getSpectrumDifference_(
      spectrum_first.begin(), spectrum_first.end(),
      spectrum_second.begin(), spectrum_second.end(),
      std::inserter(spectrum_first_diff, spectrum_first_diff.begin()));
      
    RichPeakSpectrum spectrum_second_diff;
    AScore::getSpectrumDifference_(
      spectrum_second.begin(), spectrum_second.end(),
      spectrum_first.begin(), spectrum_first.end(),
      std::inserter(spectrum_second_diff, spectrum_second_diff.begin()));
      
    site_determining_ions[0] = spectrum_first_diff;
    site_determining_ions[1] = spectrum_second_diff;
    
    site_determining_ions[0].sortByPosition();
    site_determining_ions[1].sortByPosition(); 
  }

  Size AScore::numberOfMatchedIons_(const RichPeakSpectrum & th, const RichPeakSpectrum & window, Size depth, double fragment_mass_tolerance, bool fragment_mass_tolerance_ppm) const
  {
    Size n = 0;
    
    for (Size i = 0; i < window.size() && i <= depth; ++i)
    {
      Size nearest_peak = th.findNearest(window[i].getMZ());
      
      if (nearest_peak < th.size())
      {
        double theo_mz = th[nearest_peak].getMZ();
        double error = abs(theo_mz - window[i].getMZ());

        if (fragment_mass_tolerance_ppm)
        {
          error = error / theo_mz * 1e6;
        }

        if (error < fragment_mass_tolerance)
        {
          ++n;
        }
      }
    }
    return n;
  }

  double AScore::peptideScore_(const std::vector<double> & scores) const
  {
    OPENMS_PRECONDITION(scores.size() == 10, "Scores vector must contain a score for every peak level."); 
    return (scores[0] * 0.5
            + scores[1] * 0.75
            + scores[2]
            + scores[3]
            + scores[4]
            + scores[5]
            + scores[6] * 0.75
            + scores[7] * 0.5
            + scores[8] * 0.25
            + scores[9] * 0.25)
           / 10.0;
  }

  vector<Size> AScore::getSites_(AASequence& without_phospho) const
  {
    vector<Size> tupel;
    String unmodified = without_phospho.toUnmodifiedString();
    for (Size i = 0; i < unmodified.size(); ++i)
    {
      if (unmodified[i] == 'Y' || unmodified[i] == 'T' || unmodified[i] == 'S')
      {
        tupel.push_back(i);
      }
    }
    return tupel;
  }

  vector<vector<Size> > AScore::computePermutations_(vector<Size> sites, Int n_phosphorylation_events) const
  {
    vector<vector<Size>  > permutations;
    
    if (n_phosphorylation_events == 0)
    {
      return permutations;
    }
    else if (n_phosphorylation_events == 1)
    {
      for (Size i = 0; i < sites.size(); ++i)
      {
        vector<Size> temp;
        temp.push_back(sites[i]);
        permutations.push_back(temp);
      }
      return permutations;
    }
    // All sites are phosphorylated? Return one permutation containing all sites at once.
    else if (sites.size() == (Size)n_phosphorylation_events)
    {
      permutations.push_back(sites);
      return permutations;
    }
    else
    // Generate all n_phosphorylation_events sized sets from sites
    {
      vector<Size> head;
      vector<vector<Size> > tail;
      
      // all permutations with first site selected
      head.push_back(sites[0]);
      vector<Size> tupel_left(++sites.begin(), sites.end());
      Int tail_phospho_sites = n_phosphorylation_events - 1;
      
      tail = computePermutations_(tupel_left, tail_phospho_sites);
      
      for (vector<vector<Size> >::iterator it = tail.begin(); it < tail.end(); ++it)
      {
        vector<Size> temp(head);
        temp.insert(temp.end(), it->begin(), it->end());
        permutations.push_back(temp);
      }

      // all permutations with first site not selected
      vector<vector<Size> > other_possibilities(computePermutations_(tupel_left, n_phosphorylation_events));
      permutations.insert(permutations.end(), other_possibilities.begin(), other_possibilities.end());
      return permutations;
    }
  }
  
  /// Computes number of phospho events in a sequence
  Size AScore::numberOfPhosphoEvents_(const String sequence) const 
  {
    Size cnt_phospho_events = 0;
    
    for (Size i = sequence.find("Phospho"); i != std::string::npos; i = sequence.find("Phospho", i + 7))
    {
      ++cnt_phospho_events;
    }
    
    return cnt_phospho_events;
  }
    
  /// Create variant of the peptide with all phosphorylations removed
  AASequence AScore::removePhosphositesFromSequence_(const String sequence) const 
  {
    String seq(sequence);
    seq.substitute("(Phospho)", "");
    AASequence without_phospho = AASequence::fromString(seq);
    
    return without_phospho;
  }
  
  /// Create theoretical spectra
  vector<RichPeakSpectrum> AScore::createTheoreticalSpectra_(const vector<vector<Size> > & permutations, const AASequence & seq_without_phospho) const
  {
    vector<RichPeakSpectrum> th_spectra;
    TheoreticalSpectrumGenerator spectrum_generator;
    
    th_spectra.resize(permutations.size());
    for (Size i = 0; i < permutations.size(); ++i)
    {
      AASequence seq(seq_without_phospho);
      Size permu = 0;
      
      for (Size as = 0; as < seq.size(); ++as)
      {
        if (as == permutations[i][permu])
        {
          seq.setModification(as, "Phospho");
          ++permu;
        }
        
        if (permu == permutations[i].size()) {
          break;
        }
      }

      // we mono-charge spectra
      spectrum_generator.addPeaks(th_spectra[i], seq, Residue::BIon, 1);
      spectrum_generator.addPeaks(th_spectra[i], seq, Residue::YIon, 1);
      th_spectra[i].setName(seq.toString());
    }
    return th_spectra;
  }
  
  vector<RichPeakSpectrum> AScore::createTheoreticalSpectra_(const AASequence & seq_without_phospho) const
  {
    vector<RichPeakSpectrum> th_spectra;
    TheoreticalSpectrumGenerator spectrum_generator;
    
    // we mono-charge spectra
    th_spectra.resize(1);
    spectrum_generator.addPeaks(th_spectra[0], seq_without_phospho, Residue::BIon, 1);
    spectrum_generator.addPeaks(th_spectra[0], seq_without_phospho, Residue::YIon, 1);
    th_spectra[0].setName(seq_without_phospho.toString());
    
    return th_spectra;
  }
  
  std::vector<RichPeakSpectrum> AScore::peakPickingPerWindowsInSpectrum_(RichPeakSpectrum &real_spectrum) const
  {
    vector<RichPeakSpectrum> windows_top10;
    
    double spect_lower_bound = floor(real_spectrum.front().getMZ() / 100) * 100;
    double spect_upper_bound = ceil(real_spectrum.back().getMZ() / 100) * 100;
    Size number_of_windows = static_cast<Size>(ceil((spect_upper_bound - spect_lower_bound) / 100));
    windows_top10.resize(number_of_windows);
    
    RichPeakSpectrum::Iterator it_current_peak = real_spectrum.begin();
    Size window_upper_bound(spect_lower_bound + 100);
    
    for (Size current_window = 0; current_window < number_of_windows; ++current_window)
    {
      RichPeakSpectrum real_window;
      while (((*it_current_peak).getMZ() <= window_upper_bound) && (it_current_peak < real_spectrum.end()))
      {
        real_window.push_back(*it_current_peak);
        ++it_current_peak;
      }
      
      real_window.sortByIntensity(true);
      for (Size i = 0; (i < 10) & (i < real_window.size()); ++i)
      {
        windows_top10[current_window].push_back(real_window[i]);        
      }
      
      window_upper_bound += 100;
    }
    return windows_top10;
  }
  
  std::vector<std::vector<double> > AScore::calculatePermutationPeptideScores_(vector<RichPeakSpectrum>& th_spectra, const vector<RichPeakSpectrum>& windows_top10, double fragment_mass_tolerance, bool fragment_mass_unit_ppm) const
  {
    //prepare peak depth for all windows in the actual spectrum
    vector<vector<double> > permutation_peptide_scores(th_spectra.size());
    vector<vector<double> >::iterator site_score = permutation_peptide_scores.begin();
    
    // for each phospho site assignment
    for (vector<RichPeakSpectrum>::iterator it = th_spectra.begin(); it < th_spectra.end(); ++it, ++site_score)
    {
      // the number of theoretical peaks (all b- and y-ions) correspond to the number of trials N
      Size N = it->size();
      site_score->resize(10);
      for (Size i = 1; i <= 10; ++i)
      {
        Size n = 0;
        for (Size current_win = 0; current_win < windows_top10.size(); ++current_win) // count matched ions over all 100 Da windows
        {
          n += numberOfMatchedIons_(*it, windows_top10[current_win], i, fragment_mass_tolerance, fragment_mass_unit_ppm);
        }
        double p = static_cast<double>(i) / 100.0;
        double cumulative_score = computeCumulativeScore_(N, n, p);
        
        //abs is used to avoid -0 score values
        (*site_score)[i - 1] = abs((-10.0 * log10(cumulative_score)));
      }
    }
    return permutation_peptide_scores;
  }
  
  std::multimap<double, Size> AScore::rankWeightedPermutationPeptideScores_(const vector<vector<double> > & peptide_site_scores) const
  {
    multimap<double, Size> ranking;
    for (Size i = 0; i != peptide_site_scores.size(); ++i)
    {
      double weighted_score = peptideScore_(peptide_site_scores[i]);
      ranking.insert(pair<double, Size>(weighted_score, i));
    }
    
    return ranking;
  }
} // namespace OpenMS
